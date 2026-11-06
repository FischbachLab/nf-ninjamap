#!/usr/bin/env python3
"""Assign perfect BAM alignments to strains with compact per-read state.

Read state is O(reads + distinct read/strain matches); alignment records are never
retained. Optional BAM outputs use a disk-backed alignment-offset index to preserve record order. CLI and abundance
calculations follow the original ninjaMap_parallel.py, including its BED voting
policy and minimum-singular-vote denominator. Coverage is calculated downstream,
not here. Only pysam and pandas (for the small final strain tables) are required.
"""

import argparse
import csv
import gzip
import logging
import os
import sqlite3
import tempfile
from collections import OrderedDict
from contextlib import ExitStack, closing
from time import perf_counter

import pandas as pd
import pysam


class Read:
    __slots__ = ('candidates', 'length', 'template_length', 'singular',
                 'has_voted', 'is_singular', 'output_candidates')

    def __init__(self, aln):
        # Most reads hit very few strains. Promote to a dict only for highly
        # ambiguous reads, keeping insertion fast without a dict per read.
        self.candidates = ()
        self.length = abs(aln.reference_length)
        self.template_length = abs(aln.template_length)
        self.singular = None
        self.has_voted = False
        self.is_singular = False
        self.output_candidates = ()


class Strain:
    __slots__ = ('name', 'weight', 'genome_size', 'singular_reads', 'escrow_reads',
                 'singular_votes', 'escrow_votes', 'singular_bases',
                 'escrow_bases', 'adjusted_weight', 'index')

    def __init__(self, name, weight, index):
        self.name = name
        self.index = index
        self.weight = weight
        self.genome_size = 0
        self.singular_reads = self.escrow_reads = 0
        self.singular_votes = self.escrow_votes = 0
        self.singular_bases = self.escrow_bases = 0
        self.adjusted_weight = 0


    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name


class AlignmentSpool:
    """Index BAM offsets on disk to reproduce the original output record order."""

    def __init__(self, path):
        self.db = sqlite3.connect(path)
        self.db.execute('PRAGMA cache_size=-8192')
        self.db.execute('PRAGMA temp_store=FILE')
        self.db.execute('PRAGMA journal_mode=OFF')
        self.db.execute('PRAGMA synchronous=OFF')
        self.db.execute('CREATE TABLE alignments (name TEXT, strain INTEGER, position INTEGER)')
        self.db.execute('CREATE TABLE plan (kind TEXT, name TEXT, strain INTEGER)')
        self.batch = []
        self.blocks = OrderedDict()

    def add(self, name, strain, position):
        self.batch.append((name, strain, position))
        if len(self.batch) == 4096:
            self.flush()

    def flush(self):
        self.db.executemany('INSERT INTO alignments VALUES (?, ?, ?)', self.batch)
        self.batch.clear()

    def finish(self):
        self.flush()
        self.db.execute('CREATE INDEX alignment_lookup ON alignments(name, strain)')

    def plan(self, kind, name, strain):
        self.db.execute('INSERT INTO plan VALUES (?, ?, ?)', (kind, name, strain))

    def alignment(self, source, position):
        # Read ahead within one BGZF block, keeping at most eight blocks in RAM.
        # BAM offsets preserve auxiliary tag types that a SAM round trip loses.
        block = position >> 16
        records = self.blocks.pop(block, {})
        self.blocks[block] = records
        if len(self.blocks) > 8:
            self.blocks.popitem(last=False)
        if position not in records:
            source.seek(position)
            current = position
            while current >> 16 == block:
                try:
                    records[current] = next(source)
                except StopIteration:
                    break
                current = source.tell()
        return records[position]

    def write_read(self, output, source, name, strain):
        for (position,) in self.db.execute(
                'SELECT position FROM alignments WHERE name=? AND strain=? ORDER BY rowid',
                (name, strain)):
            output.write(self.alignment(source, position))

    def close(self):
        self.db.close()


def read_name(aln):
    # Preserve the old naming rule: strip a suffix only for exactly one slash.
    name = aln.query_name
    if name.count('/') == 1:
        name = name.split('/', 1)[0]
    return name + ('__fwd' if aln.is_read1 else '__rev')


def mate_name(name):
    return name[:-3] + ('rev' if name.endswith('fwd') else 'fwd')


def is_perfect(aln):
    if aln.is_unmapped:
        return False
    return (aln.get_tag('NM') == 0 and
            aln.query_length == aln.get_overlap(aln.reference_start, aln.reference_end))


def load_metadata(path):
    bins, indices, strains = {}, {}, []
    with open(path, newline='') as handle:
        rows = csv.reader(handle)
        next(rows, None)
        for row in rows:
            name, _, weight, _, _, _, contig, length = row
            name = os.path.basename(name)
            if name not in indices:
                indices[name] = len(strains)
                strains.append(Strain(name, float(weight), len(strains)))
            index = indices[name]
            strains[index].genome_size += int(length)
            bins[contig] = index
    return bins, strains


def alignments_with_offsets(bam, track_offsets):
    if not track_offsets:
        for aln in bam.fetch(until_eof=True):
            yield None, aln
        return
    while True:
        position = bam.tell()
        try:
            aln = next(bam)
        except StopIteration:
            return
        yield position, aln


def load_reads(path, bins, threads, spool=None):
    # Nonperfect reads need only a name; the two collections remain disjoint.
    reads = {}
    imperfect = set()
    with pysam.AlignmentFile(path, 'rb', threads=threads) as bam:
        reference_bins = [bins.get(name) for name in bam.references]
        for position, aln in alignments_with_offsets(bam, spool is not None):
            name = read_name(aln)
            read = reads.get(name)
            perfect = is_perfect(aln)
            if read is None:
                if not perfect:
                    imperfect.add(name)
                    continue
                imperfect.discard(name)
                read = Read(aln)
                reads[name] = read
            else:
                # The original uses lengths from the last alignment, even if
                # that alignment is imperfect.
                read.length = abs(aln.reference_length or 0)
                read.template_length = abs(aln.template_length)
            if perfect:
                index = reference_bins[aln.reference_id]
                if index is None:
                    raise ValueError('Contig absent from bin map: ' + aln.reference_name)
                if spool is not None:
                    spool.add(name, index, position)
                if isinstance(read.candidates, dict):
                    read.candidates[index] = None
                elif index not in read.candidates:
                    read.candidates += (index,)
                    if len(read.candidates) > 4:
                        read.candidates = dict.fromkeys(read.candidates)
    for read in reads.values():
        read.candidates = tuple(read.candidates)
    return reads, len(reads) + len(imperfect)


def load_masked(path, bam_path, reads):
    masked = {}
    if not path or not os.path.isfile(path) or not os.path.getsize(path):
        return masked
    with pysam.AlignmentFile(bam_path, 'rb') as bam, open(path) as bed:
        for line in bed:
            if line.startswith('#') or not line.strip():
                continue
            contig, start, end = line.split()[:3]
            for aln in bam.fetch(contig, int(start), int(end)):
                name = read_name(aln)
                if name in reads:
                    masked[name] = None
    return masked


def intersection(read, mate, strains):
    # Match the original name-based object hashes and set iteration order. This
    # preserves both floating-point summation order and the order of vote rows.
    return tuple(strain.index for strain in
                 set(strains[i] for i in read.candidates) &
                 set(strains[i] for i in mate.candidates))


def candidates_for(read, mate, strains):
    if mate is not None:
        common = intersection(read, mate, strains)
        if common:
            return common
    return read.candidates


VOTE_HEADER = ('Read_Name', 'Strain_Name', 'cSingular_Vote', 'cEscrow_Vote',
               'Num_Strains_Aligned', 'was_discarded', 'has_voted',
               'is_singular', 'mate_has_perfect_aln')


def assign_votes(reads, strains, masked, vote_path, coverage, spool=None):
    initial_singular = 0
    for name, read in reads.items():
        read.is_singular = len(read.candidates) == 1 and name not in masked
        initial_singular += read.is_singular

    # Keep the legacy policy: a masked read skips initiating a pair, but its
    # unmasked mate may still assign that pair. Changing this needs a separate
    # biological decision rather than a performance refactor.
    escrow = {name: reads[name] for name in masked}
    singular_count = 0
    with gzip.open(vote_path, 'wt', newline='') as handle:
        writer = csv.writer(handle, lineterminator='\n')
        writer.writerow(VOTE_HEADER)
        for name, read in reads.items():
            if read.has_voted or name in masked:
                continue
            mate_key = mate_name(name)
            mate = reads.get(mate_key)
            common = (() if mate is None else
                      intersection(read, mate, strains))
            if len(common) != 1:
                escrow[name] = read
                continue
            index = common[0]
            strain = strains[index]
            strain.singular_bases += read.template_length
            strain.singular_reads += 2
            strain.singular_votes += 2
            singular_count += 2
            for key, item in ((name, read), (mate_key, mate)):
                if spool is not None and coverage:
                    spool.plan('singular', key, index)
                item.singular = index
                item.has_voted = True
                writer.writerow((key, strain.name, 1, 0, 1, False, True, item.is_singular, True))

        kept = discarded = 0
        for name, read in escrow.items():
            if read.has_voted:
                continue
            if name in masked:
                read.has_voted = True
                discarded += 1
                continue
            mate_key = mate_name(name)
            approved = candidates_for(read, escrow.get(mate_key), strains)
            weights = []
            for index in approved:
                strain = strains[index]
                strain.adjusted_weight = strain.singular_votes / strain.weight
                weights.append((index, strain.adjusted_weight))
            total_weight = 0
            for _, weight in weights:
                total_weight += weight
            to_discard = total_weight <= 0
            if to_discard:
                discarded += 1
            else:
                kept += 1
                read.has_voted = True
                if coverage:
                    read.output_candidates = tuple(approved)
            cumulative = 0
            for index, weight in weights:
                strain = strains[index]
                vote = weight / total_weight if not to_discard else 0
                if not to_discard:
                    strain.escrow_votes += vote
                    strain.escrow_reads += 1
                    strain.escrow_bases += read.length
                    cumulative += vote
                if spool is not None and coverage and not to_discard and not read.is_singular:
                    spool.plan('escrow', name, index)
                writer.writerow((name, strain.name, 0, vote, len(approved),
                                 to_discard, read.has_voted, read.is_singular,
                                 mate_key in reads))
            if not to_discard and not (cumulative < 0.001 or 0.999 < cumulative < 1.001):
                raise RuntimeError('Voter fraud detected for ' + name)
    return initial_singular, singular_count, len(escrow), kept, discarded


def write_bams(args, prefix, reads, strains, spool):
    if spool is None:
        return
    logging.info('Writing optional BAMs in legacy record order from the disk spool ...')
    with ExitStack() as stack:
        bam = stack.enter_context(pysam.AlignmentFile(args.bamfile, 'rb'))
        outputs = {}
        suffixes = []
        if args.coverage:
            suffixes += ['singular', 'escrow']
        if args.debug:
            suffixes += ['ninjaMap.false_positives', 'ninjaMap.true_positives']
        for suffix in suffixes:
            outputs[suffix] = stack.enter_context(
                pysam.AlignmentFile(prefix + '.' + suffix + '.bam', 'wb', template=bam))
        if args.debug:
            for name, read in reads.items():
                if len(read.candidates) == 1:
                    index = read.candidates[0]
                    suffix = ('ninjaMap.true_positives' if strains[index].name == args.truth
                              else 'ninjaMap.false_positives')
                    spool.write_read(outputs[suffix], bam, name, index)
        if args.coverage:
            for kind, name, index in spool.db.execute('SELECT kind, name, strain FROM plan ORDER BY rowid'):
                spool.write_read(outputs[kind], bam, name, index)


def check_votes(path, temp_parent=None):
    """Independent exact per-read audit, with bounded RAM and disk-backed totals."""
    with tempfile.TemporaryDirectory(prefix='ninjaMap-audit-', dir=temp_parent) as tmp:
        # A sqlite3 connection context only manages the transaction. Close it
        # explicitly before TemporaryDirectory removes files, including on NFS.
        with closing(sqlite3.connect(os.path.join(tmp, 'totals.sqlite'))) as db:
            db.execute('PRAGMA cache_size=-8192')
            db.execute('PRAGMA journal_mode=OFF')
            db.execute('PRAGMA synchronous=OFF')
            db.execute('CREATE TABLE totals (name TEXT PRIMARY KEY, votes REAL) WITHOUT ROWID')
            sql = ('INSERT INTO totals VALUES (?, ?) ON CONFLICT(name) DO UPDATE '
                   'SET votes=totals.votes+excluded.votes')
            with gzip.open(path, 'rt', newline='') as handle:
                batch = []
                for row in csv.DictReader(handle):
                    batch.append((row['Read_Name'], round(float(row['cSingular_Vote']) +
                                                         float(row['cEscrow_Vote']), 5)))
                    if len(batch) == 16384:
                        db.executemany(sql, batch)
                        batch.clear()
                db.executemany(sql, batch)
            with closing(db.execute('SELECT 1 FROM totals WHERE votes > 1.001 OR '
                                    '(votes > 0.001 AND votes < 0.999) LIMIT 1')) as cursor:
                return cursor.fetchone() is not None


# Preserve the original declared CSV schema explicitly, independent of pandas
# DataFrame.add/column-union behavior in the runtime environment.
STRAIN_STATS_COLUMNS = (
    'Genome_Size', 'Read_Fraction', 'Singular_Strain_Weight',
    'Total_Singular_Reads', 'Total_Singular_Votes', 'Singular_Read_Vote_Ratio',
    'Singular_Fraction_of_Singular_Reads', 'Singular_Coverage',
    'Total_Escrow_Reads', 'Total_Escrow_Votes', 'Escrow_Read_Vote_Ratio',
    'Fraction_of_all_Escrow_Reads', 'Escrowed_Cov',
)


def write_tables(prefix, strains, total, perfect, counts, threshold):
    initial_singular, singular, escrow, kept, discarded = counts
    if threshold > 0:
        for strain in strains:
            if strain.singular_votes <= threshold:
                total -= strain.singular_votes
    stats, abundance = [], []
    for strain in sorted(strains, key=lambda item: item.name):
        fraction = ((strain.escrow_votes + strain.singular_votes) * 100 / total
                    if total and strain.singular_votes > threshold else 0)
        stats.append({
            'Genome_Size': strain.genome_size,
            'Read_Fraction': fraction,
            'Singular_Strain_Weight': strain.adjusted_weight,
            'Total_Singular_Reads': strain.singular_reads,
            'Total_Singular_Votes': strain.singular_votes,
            'Singular_Read_Vote_Ratio': (strain.singular_votes / strain.singular_reads
                                        if strain.singular_reads else 0),
            'Singular_Fraction_of_Singular_Reads': (strain.singular_reads / initial_singular
                                                    if initial_singular else 0),
            'Singular_Coverage': strain.singular_bases,
            'Total_Escrow_Reads': strain.escrow_reads,
            'Total_Escrow_Votes': strain.escrow_votes,
            'Escrow_Read_Vote_Ratio': (strain.escrow_votes / strain.escrow_reads
                                      if strain.escrow_reads else 0),
            'Fraction_of_all_Escrow_Reads': strain.escrow_reads / escrow if escrow else 0,
            'Escrowed_Cov': strain.escrow_bases,
        })
        abundance.append({'Read_Fraction': fraction})
    names = [strain.name for strain in sorted(strains, key=lambda item: item.name)]
    for rows, suffix in ((stats, 'strain_stats'), (abundance, 'abundance')):
        columns = STRAIN_STATS_COLUMNS if suffix == 'strain_stats' else ('Read_Fraction',)
        frame = pd.DataFrame(rows, index=names, columns=columns)
        # The initial empty DataFrame.add promotes every legacy numeric column
        # to float, including counts and lengths, even for one strain.
        frame = frame.astype(float)
        frame.to_csv(prefix + '.ninjaMap.' + suffix + '.csv', index_label='Strain_Name')
    with open(prefix + '.ninjaMap.read_stats.csv', 'w', newline='') as handle:
        writer = csv.writer(handle, lineterminator='\n')
        writer.writerow(('File_Name', 'Reads_Aligned', 'Reads_wPerfect_Aln',
                         'Reads_wSingular_Votes', 'Reads_wEscrowed_Votes',
                         'Discarded_Reads_w_Perfect_Aln'))
        writer.writerow((prefix, total, perfect, singular, kept, discarded))


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-bam', dest='bamfile', help='Input BAM; any sort order.')
    p.add_argument('-bin', dest='binmap', help='Eight-column strain index CSV.')
    p.add_argument('-outdir')
    p.add_argument('--audit-only', metavar='VOTES.CSV.GZ',
                   help='Audit an existing votes file without reprocessing a BAM.')
    p.add_argument('--tmpdir', '-tmpdir', default=None,
                   help='SQLite scratch directory; defaults to Python temporary storage. '
                        'Use a local disk directory when outputs are on EFS/NFS.')
    p.add_argument('-bed', dest='bedfile', default='')
    p.add_argument('-prefix', default='')
    p.add_argument('-threads', type=int, default=1,
                   help='BAM decompression threads (default: 1). No process pools.')
    p.add_argument('-debug', action='store_true')
    p.add_argument('-truth', default=False)
    p.add_argument('-mbq', dest='min_base_qual', type=int, default=20,
                   help='Accepted for compatibility; coverage is calculated downstream.')
    p.add_argument('-coverage', action='store_true', help='Write singular and escrow BAMs.')
    p.add_argument('-msv', dest='min_singular_vote', type=int, default=0)
    p.add_argument('-version', action='store_true', help='Accepted for CLI compatibility.')
    return p


def main(argv=None):
    argument_parser = parser()
    args = argument_parser.parse_args(argv)
    if not args.audit_only:
        missing = [flag for flag, value in (('-bam', args.bamfile),
                                            ('-bin', args.binmap),
                                            ('-outdir', args.outdir)) if not value]
        if missing:
            argument_parser.error('the following arguments are required: ' + ', '.join(missing))
    if args.threads < 1:
        argument_parser.error('-threads must be positive')
    if args.tmpdir:
        os.makedirs(args.tmpdir, exist_ok=True)
    started = perf_counter()
    logging.basicConfig(level=logging.INFO, format='%(asctime)s\t[%(levelname)s]:\t%(message)s')
    if args.audit_only:
        logging.info('Auditing existing votes without reprocessing the BAM: %s', args.audit_only)
        if check_votes(args.audit_only, args.tmpdir):
            raise RuntimeError('Independent voter audit detected inconsistent votes')
        logging.info('Independent voter audit passed in %.2f seconds', perf_counter() - started)
        return
    os.makedirs(args.outdir, exist_ok=True)
    prefix = os.path.join(args.outdir, args.prefix or os.path.basename(args.bamfile).split('.')[0])
    logging.info('Loading strain metadata and perfect read candidates ...')
    bins, strains = load_metadata(args.binmap)
    with ExitStack() as stack:
        spool = None
        if args.coverage or args.debug:
            temporary = stack.enter_context(tempfile.TemporaryDirectory(
                prefix='ninjaMap-alignments-', dir=args.tmpdir))
            spool = AlignmentSpool(os.path.join(temporary, 'alignments.sqlite'))
            stack.callback(spool.close)
        reads, total = load_reads(args.bamfile, bins, args.threads, spool)
        if not total:
            raise SystemExit('No aligned reads found!')
        if not reads:
            raise SystemExit('No perfect aligned reads found!')
        if spool is not None:
            spool.finish()
        perfect = len(reads)
        logging.info('%d reads; %d with perfect alignments; %d strains', total, perfect, len(strains))
        masked = load_masked(args.bedfile, args.bamfile, reads)
        logging.info('Found %d masked reads. Assigning singular and escrow votes ...', len(masked))
        vote_path = prefix + '.ninjaMap.votes.csv.gz'
        counts = assign_votes(reads, strains, masked, vote_path, args.coverage, spool)
        write_bams(args, prefix, reads, strains, spool)
        del reads, masked, bins
    write_tables(prefix, strains, total, perfect, counts, args.min_singular_vote)
    logging.info('Running independent disk-backed voter audit ...')
    if check_votes(vote_path, args.tmpdir):
        raise RuntimeError('Independent voter audit detected inconsistent votes')
    logging.info('Completed in %.2f seconds; singular=%d, escrow kept=%d, discarded=%d',
                 perf_counter() - started, counts[1], counts[3], counts[4])


if __name__ == '__main__':
    main()
