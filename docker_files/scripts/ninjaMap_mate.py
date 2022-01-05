#!/usr/bin/env python3

# Each read gets 1 vote for candidate strains. 
# The more strains it maps to the lower the value of it's vote.
# Note: we are not calculating coverage, just read assignments.
# MAJOR ASSUMPTION 1: We will always have a complete genome for organisms (exact strain) that we seek in the sample.

import argparse
import gzip
import logging
# import multiprocessing
import numpy as np
import os
import pybedtools
import pysam
import pandas as pd
import re
import sys

from collections import defaultdict, Counter
from time import perf_counter as timer

start = timer()
###############################################################################
# Functions for making this more robust
###############################################################################
def human_time(time):
    time = abs(time)
    day = time // (24 * 3600)
    time = time % (24 * 3600)
    hour = time // 3600
    time %= 3600
    minutes = time // 60
    time %= 60
    seconds = time
    time_str = format('%02d:%02d:%02d:%02d'%(day,hour,minutes,seconds))
    return time_str

def intersection(list1, list2):
    return list(set(list1) & set(list2))

# def bam_is_empty(fn):
#     if os.path.getsize(fn) > 1000000:
#         return False

#     bam = pysam.Samfile(fn, check_sq=False)
#     try:
#         bam.next()
#         return False
#     except StopIteration:
#         return True

def sort_and_index(file_name, cores=4, by='Coord'):
    """ Sorts and indexes a bam file by coordinates.
    """
    if by == 'Coord':
        sorted_name = file_name.replace('.bam', '') + '.sortedByCoord.bam'
        # pysam sort multithreading support doesn't work
        # pysam.sort('-@',cores,'-o',sorted_name, file_name)
        pysam.sort('-o',sorted_name, file_name)
        
    elif by == 'Name':
        sorted_name = file_name.replace('.bam', '') + '.sortedByName.bam'
        pysam.sort('-n','-o',sorted_name, file_name)
    else:
        raise Exception("Bam file can only be sorted by 'Coord' or 'Name'.")

    pysam.index(sorted_name)
    os.remove(file_name)
    return sorted_name

###############################################################################
# Side Project
###############################################################################

# def predict_alignment():
#     feature_cols = ["align_len","query_len","aln_cov","quality","perc_id","aln_score","mate_score","mismatches","gap_open","gap_ext","is_dup","is_primary","is_supp"]
#     read the model pickle, predict and write outcome.

###############################################################################
# Setup Input and Script Usage
###############################################################################

usage = """
    USAGE: 
    python ninjaMap.py \
-bam input_bamfile \
-bin tab-delimited file with Col1= contig name and Col2=Bin/Strain name \
-out abundance table output
-log logfile.txt
    """

p = argparse.ArgumentParser(   
    formatter_class=argparse.RawTextHelpFormatter,
    add_help=True,
    usage=argparse.SUPPRESS,
    description="""Description:
This script will calculate the abundance of a strain in a defined microbial community. 
Usage: ninjaMap.py -bam name_sorted.bam -bin contig_strain_assignments.tsv -out abundance_table_output.tsv
""",
    epilog="""Examples:
python ninjaMap.py -bin contig_names_bin_map.txt -bam Bacteroides-sp-9-1-42FAA/Bacteroides-sp-9-1-42FAA.processed.sortedByCoord.bam -prefix Bacteroides-sp-9-1-42FAA    
""")
# Required
p.add_argument('-bam', dest='bamfile', action='store', type=str, required = True,
                help='name sorted bam file.')
p.add_argument('-fasta', dest='fastafile', action='store', type=str, required = True,
                help='database fasta file')
p.add_argument('-bin', dest='binmap', action='store', type=str, required = True,
                help='tab-delimited file with Col1= contig name and Col2=Bin/Strain name')
p.add_argument('-outdir', dest='outdir', action='store', type=str, required = True,
                help='output directory')
# Optional
p.add_argument('-prefix', dest='prefix', action='store', type=str,
                help='output prefix')
p.add_argument('-cores', dest='cores', action='store', default=8,
                help='number of cores available for use')
p.add_argument('-test', dest='test', action='store_true', default=False,    
                help='save intermediate false positives bam file')
p.add_argument('-truth', dest='truth', action='store', default=False,    
                help='If using test, please provide one strain name that you would like to track.')
p.add_argument('-mbq', dest='min_base_qual', action='store', default=20, type=int,    
                help='minimum read base quality to consider for coverage calculations.')

args = vars(p.parse_args())
# test
# bamfile_name = "/Users/sunit.jain/Research/SyntheticCommunities/ReadAlignment/Testing/Mismaps/Bacteroides-coprophilus-DSM-18228/Bacteroides-coprophilus-DSM-18228.processed.bam"
# abundance_output_file = 'B_coprophilius.ninjaMap.v1.abundance.tsv'
bamfile_name = args['bamfile']
fastafile_name = args['fastafile']
binmap_file = args['binmap']
output_dir = args['outdir']

cores = args['cores']
os.environ['NUMEXPR_MAX_THREADS'] = f'{cores}'
os.makedirs(output_dir, exist_ok=True)
os.makedirs(f'{output_dir}/tmp', exist_ok=True)

if not args['prefix']:
    default_prefix = os.path.basename(bamfile_name).split('.')[0]
    prefix = os.path.join(output_dir, default_prefix)
else:
    prefix = os.path.join(output_dir, args['prefix'])

abundance_output_file = prefix +'.ninjaMap.abundance.csv'
stats_file = prefix +'.ninjaMap.read_stats.csv'
vote_file = prefix +'.ninjaMap.votes.csv.gz'
fraud_file = prefix +'.ninjaMap.fraud_votes.txt'
strain_stats_file = prefix +'.ninjaMap.strain_stats.csv'
logfile = prefix +'.ninjaMap.log.txt'

if args['min_base_qual']:
    min_base_qual = args['min_base_qual']
else:
    min_base_qual = 20

logging.basicConfig(
    # filename=logfile, 
    # filemode='w+', 
    level=logging.DEBUG,
    format='%(asctime)s\t[%(levelname)s]:\t%(message)s')

logging.info('Started')

singularBamfile = prefix +'.singularAln.ninjaMap.bam'
escrowBamfile = prefix +'.escrowAln.ninjaMap.bam'
tmp_bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
singularBam = pysam.AlignmentFile(singularBamfile, "wb", template=tmp_bamfile)
escrowBam = pysam.AlignmentFile(escrowBamfile, "wb", template=tmp_bamfile)

if args['test']:
    if (args['truth']):
        true_strain = args['truth']
        false_positives = pysam.AlignmentFile(prefix + '.ninjaMap.false_positives.bam', "wb", template=tmp_bamfile)
        true_positives = pysam.AlignmentFile(prefix + '.ninjaMap.true_positives.bam', "wb", template=tmp_bamfile)
    else:
        logging.critical('True strain name required for testing. None provided. Exiting.')
        sys.exit(1)

tmp_bamfile.close()

logging.info(f'Exclusive read alignments will be written to {singularBamfile}')
logging.info(f'Shared read alignments will be written to {escrowBamfile}')
###############################################################################
# Classes
###############################################################################
class Strains:
    total_genome_size = 0
    total_strains = 0
    total_uniquely_covered_bases = 0
    net_promiscuity = 0
    
    def __init__(self, strain_name):
        self.name = strain_name
        self.num_singular_reads = 0
        self.num_escrow_reads = 0
        self.num_contigs = 0
        self.genome_size = 0
        self.cum_primary_votes = 0
        self.cum_escrow_votes = 0
        self.uniqueness_score = 0
        self.percent_coverage = 0
        self.depth_variance = 0
        
        self.breadth_of_coverage = 0
        self.depth_of_coverage = 0

        self.uniquely_covered_bases = set()
        self.num_uniquely_covered_bases = 0
        self.uniquely_covered_depth = 0
        self.adj_primary_wt = 0
        self.singular_depth_var = 0
        
        self.escrow_covered_bases = set()
        self.num_escrow_covered_bases = 0
        self.escrow_covered_depth = 0
        self.escrow_depth_var = 0
        
        self.aln_norm_abundance = 0
        self.genome_norm_abundance = 0
        self.adjusted_votes = 0
        self.escrow_vote_conversion_rate = 0
        self.singular_vote_conversion_rate = 0
        
        self.contigs = defaultdict(int)
        self.singular_bin = defaultdict(int)
        self.escrow_bin = defaultdict(int)

    def __hash__(self):
        return hash(str(self.name))

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)
        
    def add_contig(self, contig_name, contig_length):
        self.contigs[contig_name] = contig_length
        self.genome_size += contig_length
        self.num_contigs += 1
        
    def add_paired_singular_vote(self, read_name, mate_name, read_template_len, read_vote = 1):
#         print(str(self) +'\t:\t' + self.name +'\t:\t' + read.unique_name +'\t:\t'+ str(read_vote))
        # self.uniquely_covered_bases += read_template_len
        # for R1
        self.singular_bin[read_name] += read_vote
        self.cum_primary_votes += read_vote 
        self.num_singular_reads += 1

        # for R2
        self.singular_bin[mate_name] += read_vote
        self.cum_primary_votes += read_vote
        self.num_singular_reads += 1
        
    
    def add_escrow_vote(self, read_name, read_template_len, read_vote):
        if read_name in self.singular_bin.keys():
            return

        self.escrow_bin[read_name] += read_vote
        self.cum_escrow_votes += read_vote
        self.num_escrow_reads += 1
        # self.escrow_covered_bases += read_template_len

    def calculate_genome_coverage(self, df, singular=True):
        if self.num_singular_reads == 0:
            return (0,0,0,0)

        coverage, depth, depth_sd, depth_variance = (0,0,0,0)

        # Scaling factor for escrow depth values
        scale = 1
        if not singular:
            if Strains.net_promiscuity > 0:
                scale = self.adj_primary_wt / Strains.net_promiscuity
            self.escrow_covered_bases = set()
            self.escrow_depth = 0
            self.num_escrow_covered_bases = 0
            self.escrow_depth_var = 0

        # filter by contig names for this strain
        contigs_df = df[df.contig.isin(self.contigs.keys())].copy()
        num_bases = contigs_df.shape[0]
        if num_bases > 0:
            contigs_df['unique_name'] = contigs_df['contig'] + '_' + contigs_df['pos'].astype(str)
            # contigs_df['unique_name'] = contigs_df[['contig', 'pos']].apply(lambda x: '_'.join(str(x)), axis = 1)
            # print(f'{contigs_df.shape}')
            # coverage = contigs_df.num_bases[contigs_df.depth > 0].sum() * 100 / self.genome_size
            # wt_depth = (contigs_df.depth * contigs_df.num_bases).sum() / self.genome_size
            # depth_variance = (contigs_df.depth * contigs_df.num_bases).var() * 100 / self.genome_size
            
            # Breadth of coverage > 0.
            coverage_breadth = contigs_df.pos[contigs_df.depth > 0].count()
            coverage = coverage_breadth * 100 / self.genome_size
            depth = contigs_df.depth.sum() / self.genome_size
            # Pad depth_array with 0 to make it equal to the genome size before calculationf SD and Coeff of Var.
            depth_array = np.concatenate([contigs_df.depth.to_numpy(),np.zeros(self.genome_size - coverage_breadth)]) * scale
            depth_sd = np.std(depth_array)
            depth_variance = depth_sd/np.mean(depth_array) # coeff of variation
            bases_covered = set(contigs_df.unique_name[contigs_df.depth > 0])
            if singular:
                self.uniquely_covered_bases = bases_covered
                self.singular_depth = depth
                self.num_uniquely_covered_bases = coverage_breadth
                self.singular_depth_var = depth_variance
                Strains.total_uniquely_covered_bases += coverage_breadth
            else:
                self.escrow_covered_bases = bases_covered
                self.escrow_depth = depth
                self.num_escrow_covered_bases = coverage_breadth
                self.escrow_depth_var = depth_variance
        
        # Allow this temporary dataframe to be garbage collected.
        del contigs_df

        return (coverage, depth, depth_sd, depth_variance)

    def _normalize_votes(self):
        # Read to Vote conversion ratio
        if self.num_singular_reads > 0:
            self.singular_vote_conversion_rate = self.cum_primary_votes/self.num_singular_reads
        else:
            self.singular_vote_conversion_rate = 0
            
        if self.num_escrow_reads > 0:
            self.escrow_vote_conversion_rate = self.cum_escrow_votes/self.num_escrow_reads
        else:
            self.escrow_vote_conversion_rate = 0

    def compile_general_stats(self):
        '''
        for each object of this class, return a pandas data frame with 
        strains as rows and number of singular and escrow votes
        '''
        if (self.num_uniquely_covered_bases == 0) and (self.num_escrow_covered_bases == 0):
            return None
        
        # Depth of cov should be added since escrow depth is already adjusted by escrow read weights.
        self.depth_of_coverage = self.escrow_depth + self.singular_depth
        
        # Breadth should be the union of Escrow and Singular coverage
        self.breadth_of_coverage = len(self.uniquely_covered_bases.union(self.escrow_covered_bases))
        self.percent_coverage = self.breadth_of_coverage * 100 / self.genome_size

        # Right now, just prioritizing reads with exclusive matches to the genome and their variance. 
        # Might not be the best way forward for cases where the exact genome is missing but multiple similar strains are present.
        self.depth_variance = self.singular_depth_var 
        self._normalize_votes()

        frac_escrow_reads = 0
        if Reads.total_escrow_reads > 0:
            frac_escrow_reads = self.num_escrow_reads/Reads.total_escrow_reads

        frac_singular_reads = 0
        if Reads.total_singular_reads_after_recruitment > 0:
            frac_singular_reads = self.num_singular_reads/Reads.total_singular_reads_after_recruitment
        
        # Dataframe
        return pd.DataFrame(
            index = [self.name],
            data  = {
                'Genome_Size' : self.genome_size,
                'Percent_Coverage' : self.percent_coverage,
                'Total_Bases_Covered' : self.breadth_of_coverage,
                'Coverage_Depth' : self.depth_of_coverage,
                'Depth_Variation' : self.depth_variance,
                'Read_Fraction' : self.read_fraction,
                'Singular_Strain_Weight' : self.adj_primary_wt,
                'Total_Singular_Reads' : self.num_singular_reads,
                'Total_Singular_Votes' : self.cum_primary_votes,
                'Singular_Read_Vote_Ratio' : self.singular_vote_conversion_rate,
                'Singular_Fraction_of_Singular_Reads' : frac_singular_reads,
                'Singular_Coverage' : self.num_uniquely_covered_bases,
                'Singular_Depth' : self.singular_depth,
                'Total_Escrow_Reads' : self.num_escrow_reads,
                'Total_Escrow_Votes' : self.cum_escrow_votes,
                'Escrow_Read_Vote_Ratio' : self.escrow_vote_conversion_rate,
                'Fraction_of_all_Escrow_Reads' : frac_escrow_reads,
                'Escrowed_Cov' : self.num_escrow_covered_bases,
                'Escrowed_Depth' : self.escrow_depth
            }
        )

    def compile_by_abundance(self):
        '''
        return a pandas data frame with 1 row x 4 columns. 
            Strain_Name,Read_Fraction, Percent_Coverage, Coverage_Depth, Depth_Coeff_Variation
        '''
        if self.read_fraction == 0:
            return None

        return pd.DataFrame(
                    index = [self.name],
                    data  = {
                        'Read_Fraction' : self.read_fraction,
                        'Percent_Coverage' : self.percent_coverage,
                        'Coverage_Depth' : self.depth_of_coverage,
                        'Depth_Coeff_Variation' : self.depth_variance
                        }
                    )
    def strain_promiscuity_adjustment(self):
        # self.adj_primary_wt = calculate_sunits_original_adjustment(self)
        self.adj_primary_wt = self.beta_adjustment()
        Strains.net_promiscuity += self.adj_primary_wt
        return self.adj_primary_wt

    # The OG
    def calculate_sunits_original_adjustment(self):
        return self.cum_primary_votes / Reads.total_reads_aligned

    def beta_adjustment(self):
        if self.num_uniquely_covered_bases == 0:
            return 0

        # 2019-11-07 Sunit's interpretation 2 of mike drop: This should be the rate of aggregation of reads for 1 strain compared to the others
        # Borrowed from the beta measure of how well a stock does compared to the rest of the sector.
        strains_reads_per_base = self.num_singular_reads / self.num_uniquely_covered_bases / self.genome_size
        
        if (Reads.total_singular_reads_after_recruitment == self.num_singular_reads):
            # all singular reads assigned were assigned to this strain.
            others_singular_reads_aligned = self.num_singular_reads
            others_uniquely_covered_bases = self.num_uniquely_covered_bases
        elif(Reads.total_singular_reads_after_recruitment > self.num_singular_reads):
            others_singular_reads_aligned = Reads.total_singular_reads_after_recruitment - self.num_singular_reads
            others_uniquely_covered_bases = Strains.total_uniquely_covered_bases - self.num_uniquely_covered_bases
        else:
            sys.exit(f"""
            {self.name} others_singular_reads_aligned can't be negative!
            Total singular reads:{Reads.total_singular_reads_after_recruitment}
            Total singular bases:{Strains.total_uniquely_covered_bases}
            Total bases:{Strains.total_genome_size}

            Strain singular reads:{self.num_singular_reads}
            Strain singular bases:{self.num_uniquely_covered_bases}
            Strain Genome Size:{self.genome_size}
            """)
        
        others_total_genomes_size = Strains.total_genome_size - self.genome_size

        others_reads_per_base = others_singular_reads_aligned / others_uniquely_covered_bases / others_total_genomes_size

        adj_primary_wt = strains_reads_per_base / others_reads_per_base
        # Case 1: self.adj_primary_wt = 1
        #       - rate of read recruitment for this strain is equal to the rate for the rest of the genomes.
        #       - No biases.
        # Case 2: self.adj_primary_wt > 1
        #       - rate of recruitment for this strain is greater than the rate for the rest of the genomes.
        #       - reads preferentially map to this strain over others.
        # Case 3: self.adj_primary_wt < 1
        #       - rate of recruitment for this strain is less than the rate for the rest of the genomes.
        #       - reads preferentially map to other strains over this strain.
        if adj_primary_wt < 0:
            sys.exit(f"""
            {self.name} adj_primary_wt can't be negative!
            Total singular reads:{Reads.total_singular_reads}
            Total singular bases:{Strains.total_uniquely_covered_bases}
            Total bases:{Strains.total_genome_size}
            
            Strain singular reads:{self.num_singular_reads}
            Strain singular bases:{self.num_uniquely_covered_bases}
            Strain Genome Size:{self.genome_size}
            """)
        return adj_primary_wt

    def calculate_read_fraction(self):
        if (self.num_uniquely_covered_bases == 0) and (self.num_escrow_covered_bases == 0):
            self.read_fraction = 0
        else:
            self.read_fraction = (self.cum_escrow_votes + self.cum_primary_votes) * 100 / Reads.total_reads_aligned
        return self.read_fraction

    @staticmethod
    def calculate_coverage(bamfile_name):
        pybedtools.set_tempdir(f'{output_dir}/tmp')
        a = pybedtools.BedTool(bamfile_name)
        df = a.genome_coverage(dz = True).to_dataframe(names=['contig','pos', 'depth'])
        pybedtools.cleanup()
        return df

class Reads:
    total_reads_aligned = 0
    reads_w_perfect_alignments = 0
    total_singular_reads = 0
    total_escrow_reads_kept = 0
    total_escrow_reads_discarded = 0
    total_escrow_reads = 0
    total_singular_reads_after_recruitment = 0

    def __init__(self, name, mate_name, read_length, template_length):
        self.name = name
        self.unique_name = name
        self.mates_unique_name = mate_name
        self.read_length = abs(read_length)
        self.template_length = abs(template_length)

        self.cum_vote = 0
        self.has_voted = False
        self.in_singular_bin = False
        self.mate_has_perfect_match = False

        self.mapped_strains = defaultdict()

    def __hash__(self):
        return hash(str(self.unique_name))

    def __eq__(self, other):
        return self.unique_name == other.unique_name

    def __ne__(self, other):
        # Not strictly necessary, but to avoid having both x==y and x!=y
        # True at the same time
        return not(self == other)

    def add_exact_match(self, strain):
        self.mapped_strains[strain] = strain.name
    
    def put_pair_in_singular_bin(self, mate):
        self.in_singular_bin = True
        mate.in_singular_bin = True

    def add_vote(self, vote_value):
        # vote_value = round(vote_value, 7)
        self.cum_vote += vote_value
        self.has_voted = True

    def is_fraud(self):
        '''
        for each object of this class, return True if cumulative votes > 1
        '''
        fraud = True
        # approx 1
        if (self.cum_vote < 1.001) and (self.cum_vote > 0.999):
            fraud = False

        # approx 0
        if (self.cum_vote < 0.001):
            fraud = False

        return fraud

    def get_voting_details(self, approved_strain_list):
        '''
        returns a list of 3 element lists, each containing: strain_name, singular vote value and escrow vote value
        '''
        vote_list = list()
        for strain in approved_strain_list:
            strain_name = ''
            escrow_votes = 0
            singular_votes = 0
            cumulative_vote = 0
            
            if self.cum_vote is not None:
                cumulative_vote = self.cum_vote

            if strain.name is not None:
                strain_name = strain.name

            if self.unique_name in strain.escrow_bin.keys():
                escrow_votes = strain.escrow_bin[self.unique_name]

            if self.unique_name in strain.singular_bin.keys():
                singular_votes = strain.singular_bin[self.unique_name]

            vote_list.append([strain_name, singular_votes, escrow_votes, cumulative_vote])
        return vote_list
    
    @staticmethod
    def choose_primary_candidate(read, mate):
        # if (read.template_length == 0) or (mate.template_length == 0):
        #     return None

        if read.mate_has_perfect_match or mate.mate_has_perfect_match :
            common_strains_list = intersection(read.mapped_strains.keys(), mate.mapped_strains.keys())
            if len(common_strains_list) == 1:
                return common_strains_list[0].name
            else:
                return None

    @staticmethod
    def is_perfect_alignment(aln):
        edit_dist = aln.get_tag('NM')
        query_len = aln.query_length
        ref_start = aln.reference_start
        ref_end = aln.reference_end
        
        # https://www.biostars.org/p/106126/
        return ((edit_dist == 0) and (query_len == aln.get_overlap(ref_start, ref_end)))
    
    @staticmethod
    def parse_read_name(aln):
        '''
        Accept: AlignmentFile object from PySam
        if read name has a '/', this is the old format. 
        strip the content after the '/', return remaining
        else, return it as is.
        '''
        try:
            key, value = aln.query_name.split("/")
        except ValueError:
            return str(aln.query_name)
        else:
            return str(key)
        
    @staticmethod
    def get_unique_read_name(aln):
        orientation = ''
        if aln.is_read1:
            orientation =  'fwd'
        else:
            orientation =  'rev'
            
        return Reads.parse_read_name(aln) +'__'+ orientation
    
    @staticmethod
    def get_unique_mate_name(aln):
        orientation = ''
        if aln.is_read1:
            orientation =  'rev'
        else:
            orientation =  'fwd'
            
        return Reads.parse_read_name(aln) +'__'+ orientation
    
    @staticmethod
    def extract_read_info(aln):
        read_name = Reads.get_unique_read_name(aln)
        mate_name = Reads.get_unique_mate_name(aln)
        read_length = aln.reference_length
        template_length = aln.template_length

        return (read_name, mate_name, read_length, template_length)
    
    @staticmethod
    def check_for_fraud(votes_file, discard_pile=1):
        fraud_status = True
        df = pd.read_csv(votes_file, index_col='Read_Name', usecols=['Read_Name', 'cSingular_Vote', 'cEscrow_Vote'])
        df['Total_Votes'] = round((df['cSingular_Vote'] + df['cEscrow_Vote']), 5)
        df = df.drop(['cSingular_Vote',  'cEscrow_Vote'], axis = 1)
        votes_df = df.groupby('Read_Name').sum()
        gt_row, gt_col = votes_df[(votes_df.Total_Votes > 1.001)].shape
        mid_row, mid_col = votes_df[(votes_df.Total_Votes > 0.001) & (votes_df.Total_Votes < 0.999)].shape
        if gt_row == 0 and mid_row == 0:
            fraud_status = False

        return fraud_status

###############################################################################
# Parse the Contig to Bin/Strain name map file.
###############################################################################

logging.info('Processing the Bin Map file: %s ...', binmap_file)
# Header = ('Strain_Name,Contig_Name,Contig_Length\n')
all_strains = defaultdict(list)
with open(binmap_file, "r") as binmap:
    for line in binmap:
        line=line.rstrip()
        contig_name, strain_name = line.split('\t')
        # all_strains[strain_name].append(line)
        all_strains[strain_name].append(contig_name)

all_strain_obj = defaultdict(int)
strains_list = sorted(list(all_strains.keys()), key=str.lower)
fasta = pysam.FastaFile(fastafile_name)
bins = defaultdict()
for strain_name in strains_list:
    strain = Strains(strain_name)
    all_strain_obj[strain_name] = strain
    for contig_name in all_strains[strain_name]:
    # for line in all_strains[strain_name]:
        # binmap_strain_name, strain.uniqueness_score, contig_name, contig_length = line.split('\t')
        strain.add_contig(contig_name, fasta.get_reference_length(contig_name))
        bins[contig_name] = strain_name
        Strains.total_genome_size += fasta.get_reference_length(contig_name)

Strains.total_strains = len(all_strains.keys())
logging.info('\t%d contigs assigned to %d strains', len(bins.keys()), Strains.total_strains)
fasta.close()

del all_strains
###############################################################################
# Parse the BAM file
###############################################################################
logging.info('Processing the BAM file: %s ...', bamfile_name)

total_reads = set()
perfect_alignment = defaultdict(lambda: defaultdict(list))
# if args['test']:
#     perfect_alignment = defaultdict(lambda: defaultdict(list))
# else:
#     perfect_alignment = defaultdict(lambda: defaultdict(int))
read_info = defaultdict()

# Read the BAM file
bamfile = pysam.AlignmentFile(bamfile_name, mode = 'rb')
for aln in bamfile.fetch(until_eof=True):
    read_name, mates_name, read_length, template_length = Reads.extract_read_info(aln)
    read_info[read_name] = (mates_name, read_length, template_length)

    if Reads.is_perfect_alignment(aln):
        strain_name = bins[aln.reference_name] # bins[contig_name] --> strain name
        perfect_alignment[read_name][strain_name].append(aln)
    
    total_reads.add(read_name)

bamfile.close()

Reads.total_reads_aligned = len(total_reads)
if Reads.total_reads_aligned == 0:
    logging.critical(f'0 reads aligned to the reference database. Please check the BAM file or the database used.',)
    sys.exit(1)

Reads.reads_w_perfect_alignments = len(perfect_alignment.keys())
logging.info('\tUsed %d reads with perfect alignments, out of %d (%7.3f%%).', 
    Reads.reads_w_perfect_alignments,
    Reads.total_reads_aligned,
    Reads.reads_w_perfect_alignments*100/Reads.total_reads_aligned
    )
del total_reads
###############################################################################
# Separate the Primary from the Escrow alignments
# Calculate the Strain abundance distribution based on the primary alignments.
###############################################################################
read_objects = defaultdict()
logging.info('Separating the Primary from the Escrow alignments ...')
for read_name in perfect_alignment.keys():
    mate_name, read_length, template_length = read_info[read_name]
    read = Reads(read_name, mate_name, read_length, template_length)
    read_objects[read_name] = read

    if mate_name in perfect_alignment.keys():
        # This means the mate had a perfect match too.
        read.mate_has_perfect_match = True

    read.num_strains = len(perfect_alignment[read_name].keys())

    for strain_name in perfect_alignment[read_name].keys():
        read.add_exact_match(all_strain_obj[strain_name])
        if args['test']:
            if read.num_strains == 1 and strain_name != true_strain:
                for aln in perfect_alignment[read_name][strain_name]:
                    # write to Bam file
                    false_positives.write(aln)
            elif read.num_strains == 1 and strain_name == true_strain:
                for aln in perfect_alignment[read_name][strain_name]:
                    # write to Bam file
                    true_positives.write(aln)

    if read.num_strains == 1:
        # Singular
        # True hit. This read gets 1 whole vote to assign to a Strain.
        read.in_singular_bin = True
        Reads.total_singular_reads += 1

if args['test']:
    false_positives.close()
    true_positives.close()

# del perfect_alignment
del read_info

logging.info('\tUsed %d reads for seed primary distribution, out of %d (%7.3f%%) reads with perfect alignments or %7.3f%% of total.', 
    Reads.total_singular_reads,
    Reads.reads_w_perfect_alignments,
    Reads.total_singular_reads*100/Reads.reads_w_perfect_alignments,
    Reads.total_singular_reads*100/Reads.total_reads_aligned
    )

if len(read_objects.keys()) != Reads.reads_w_perfect_alignments:
    logging.critical('Read %d reads with perfect alignments, but created %d read objects. There is something fishy going on here...',
    Reads.reads_w_perfect_alignments,
    len(read_objects.keys())
    )    
    sys.exit(1)

###############################################################################
# Processing the reads that mapped exclusively to a single strain
###############################################################################
votes = gzip.open(vote_file, 'wt')
fraud = open(fraud_file, 'wt')
# Header
votes.write("Read_Name,Strain_Name,cSingular_Vote,cEscrow_Vote,Num_Strains_Aligned,was_discarded,has_voted,is_singular,mate_has_perfect_aln\n")
escrow_read_objects = defaultdict()
singular_fraud_alert = False
num_singular_fraud_reads = 0
for name, read in read_objects.items():
    if read.has_voted:
        continue
    
    # if read.mates_unique_name in read_objects.keys():
    if read.mate_has_perfect_match:
        # Means both reads in a mate are prefect alignments to some strain(s); may not be the same
        mate = read_objects[read.mates_unique_name]
        strain_name = Reads.choose_primary_candidate(read, mate)

        if strain_name is not None:
            strain = all_strain_obj[strain_name]
            strain.add_paired_singular_vote(read.unique_name, mate.unique_name, read.template_length, 1)
            read.add_vote(1)
            mate.add_vote(1)

            # Write to Singular BAM file
            for aln in perfect_alignment[read.unique_name][strain.name]:
                singularBam.write(aln)
            for aln in perfect_alignment[mate.unique_name][strain.name]:
                singularBam.write(aln)

            # Reads.total_singular_reads_in_pairs += 2
            Reads.total_singular_reads_after_recruitment += 2

            votes.write(read.name + ',' + strain.name + ',' + str(strain.singular_bin[read.unique_name]) + ',' + 
                str(strain.escrow_bin[read.unique_name]) + ',' + '1,False,'+ 
                str(read.has_voted)+','+str(read.in_singular_bin)+','+str(read.mate_has_perfect_match)+'\n')
            
            votes.write(mate.name + ',' + strain.name + ',' + str(strain.singular_bin[mate.unique_name]) + ',' + 
                str(strain.escrow_bin[mate.unique_name]) + ',' + '1,False,'+ 
                str(mate.has_voted)+','+str(mate.in_singular_bin)+','+str(read.mate_has_perfect_match)+'\n')
        else:
            escrow_read_objects[name] = read
    else:
        escrow_read_objects[name] = read
    
    if read.is_fraud():
        singular_fraud_alert = True
        fraud.write(read.name+'\t'+str(read.cum_vote)+'\n')
        num_singular_fraud_reads += 1

singularBam.close()

Reads.total_escrow_reads = len(escrow_read_objects.keys())
logging.info('\t%d reads will be used for singular alignment strain abundance out of %d (%7.3f%%) reads with perfect alignments (and recruited mates) or %7.3f%% of total aligned.',
    Reads.total_singular_reads_after_recruitment,
    Reads.reads_w_perfect_alignments,
    Reads.total_singular_reads_after_recruitment*100/Reads.reads_w_perfect_alignments,
    Reads.total_singular_reads_after_recruitment*100/Reads.total_reads_aligned)
logging.info('\t%d reads for escrow alignment strain abundance out of %d (%7.3f%%) reads with perfect alignments or %7.3f%% of total aligned.',
    Reads.total_escrow_reads,
    Reads.reads_w_perfect_alignments,
    Reads.total_escrow_reads * 100 / Reads.reads_w_perfect_alignments,
    Reads.total_escrow_reads * 100 / Reads.total_reads_aligned)

if Reads.reads_w_perfect_alignments != (Reads.total_escrow_reads + Reads.total_singular_reads_after_recruitment):
    logging.critical('Read %d read objects with perfect alignments, but created %d read objects as singular and escrow total (%d + %d). There is something fishy going on here...',
    Reads.reads_w_perfect_alignments,
    (Reads.total_escrow_reads + Reads.total_singular_reads_after_recruitment),
    Reads.total_singular_reads_after_recruitment,
    Reads.total_escrow_reads)
    sys.exit(1)

del read_objects
###############################################################################
# Calculate unique number of bases covered for each genome
###############################################################################
logging.info('Computing depth and coverage for each strain in the database based on singular alignments ...')
# i = 0
# for name, strain in all_strain_obj.items():
#     i += 1
#     logging.info("\t[%d/%d] Searching for exclusive support for :\t%s",i, Strains.total_strains, name)
#     my_cov = strain.calculate_singular_coverage(bamfile_name, fastafile_name)
#     logging.info(f"\t... {my_cov*100}x\n")

# Use bedtools to calculate genomecoverage from the given bam file

# pysam.sort('-o','tmp.bam', singularBamfile)
# os.rename('tmp.bam',singularBamfile)
# pysam.index(singularBamfile)
sorted_singularBamfile = sort_and_index(singularBamfile)
singular_bed_coverage = Strains.calculate_coverage(sorted_singularBamfile)

## BEGIN: Partially successful attempt at multiprocessing ##
# def calculate_genome_coverage(strain_obj):
#     (coverage, wt_depth, depth_sd, depth_variance) = strain_obj.calculate_genome_coverage(singular_bed_coverage)
#     if coverage > 0:
#         print(f"{strain_obj.name}: Singular Coverage = {coverage}%; Weighted Depth = {wt_depth}x; Depth StDev = {depth_sd}; Variation in Depth = {depth_variance}%")
#         return (strain_obj.name, coverage, wt_depth, depth_sd, depth_variance)
## END: Partially successful attempt at multiprocessing ##

i = 0
singular_depth_file = prefix +'.ninjaMap.singular_depth.csv'
with open(singular_depth_file, 'w+') as sd_handle:
    sd_handle.write(f"strain\tcoverage\tdepth\tdepth_sd\tdepth_coeff_variance\n")
    ## BEGIN: Partially successful attempt at multiprocessing ##
    # genome_coverage = []
    # with multiprocessing.Pool(7) as p:
    #     genome_coverage.append(p.map(calculate_genome_coverage, all_strain_obj.values()))

    # print(f'{genome_coverage}')

    # for output in genome_coverage[0]:
    #     if output is not None:
    #         (name, coverage, wt_depth, depth_sd, depth_variance) = output
    #         logging.info(f"{name}: Singular Coverage = {coverage}%; Weighted Depth = {wt_depth}x; Depth StDev = {depth_sd}; Variation in Depth = {depth_variance}%")
    #         sd_handle.write(f"{name}\t{coverage}\t{wt_depth}\t{depth_sd}\t{depth_variance}\n")
    ## END: Partially successful attempt at multiprocessing ##
    for name, strain in all_strain_obj.items():
        i += 1
        logging.info("\t[%d/%d] Searching for exclusive support for :\t%s",i, Strains.total_strains, name)
        (coverage, wt_depth, depth_sd, depth_variance) = strain.calculate_genome_coverage(singular_bed_coverage)
        if coverage > 0:
            logging.info(f"\t... Calculated Singular Coverage = {coverage}%; Weighted Depth = {wt_depth}x; Depth StDev = {depth_sd}; Variation in Depth = {depth_variance}%")
            sd_handle.write(f"{strain.name}\t{coverage}\t{wt_depth}\t{depth_sd}\t{depth_variance}\n")

###############################################################################
# Use the strain abundance distribution based on Singular alignments to weight
# the Escrow abundance assignments.
###############################################################################
logging.info('Assigning escrow reads based on singular alignment strain abundance ...')
logging.info('\t and making sure there is no voter fraud ...')

escrow_fraud_alert = False
num_fraud_reads = 0
if Reads.total_escrow_reads > 0:
    for name, read in escrow_read_objects.items():
        read_vote_value = 0
        total_primary_wts = 0
        to_discard = False
        common_strains_list = list()

        # If both reads in a pair have perfect alignments,
        # shortlist the escrow vote spread to the intersection of 
        # the list of strains that each pair aligns to.
        if read.mates_unique_name in escrow_read_objects.keys():
            mate = escrow_read_objects[read.mates_unique_name]
            common_strains_list = intersection(read.mapped_strains.keys(), mate.mapped_strains.keys())
            if len(common_strains_list) == 0:
                common_strains_list = read.mapped_strains.keys()
        else:
            common_strains_list = read.mapped_strains.keys()

        if not read.has_voted:
            # Allow voting ONLY if the read it hasn't voted already
            # discard = read.calculate_escrow_abundance()
            # for strain in read.mapped_strains.keys():

            for strain in common_strains_list:
                if strain.adj_primary_wt:
                    total_primary_wts += strain.adj_primary_wt
                else:
                    total_primary_wts += strain.strain_promiscuity_adjustment()

            if total_primary_wts > 0:
                Reads.total_escrow_reads_kept += 1
                # Spread 1 vote/read based on the singular weight of strains aligned to
                for strain in common_strains_list:
                    read_vote_value = 1 * strain.adj_primary_wt / total_primary_wts
                    strain.add_escrow_vote(read.unique_name, read.read_length, read_vote_value)
                    read.add_vote(read_vote_value)
                    # Write to Escrow BAM file
                    for aln in perfect_alignment[read.unique_name][strain.name]:
                        escrowBam.write(aln)
            else:
                # not enough primary votes
                Reads.total_escrow_reads_discarded += 1
                to_discard = True

        voting_details = read.get_voting_details(common_strains_list)
        for voting_detail_sublist in voting_details:
            strain_name, singular_vote_count, escrow_vote_count, cumulative_vote = voting_detail_sublist
            votes.write(read.name + ',' + strain_name + ',' + str(singular_vote_count) + ',' + 
                        str(escrow_vote_count) + ',' + str(len(common_strains_list)) + ',' + str(to_discard) +','+ 
                        str(read.has_voted)+','+str(read.in_singular_bin)+','+ str(read.mate_has_perfect_match)+'\n')
        if read.is_fraud():
            escrow_fraud_alert = True
            fraud.write(read.name+'\t'+str(read.cum_vote)+'\n')
            num_fraud_reads += 1
votes.close()
fraud.close()
escrowBam.close()

del perfect_alignment
del escrow_read_objects

logging.info('\tUsed %d reads for escrow distribution, out of %d (%7.3f%%) reads with perfect alignments or %7.3f%% of total.',
            Reads.total_escrow_reads,
            Reads.reads_w_perfect_alignments,
            Reads.total_escrow_reads*100/Reads.reads_w_perfect_alignments,
            Reads.total_escrow_reads*100/Reads.total_reads_aligned
    )
if Reads.total_escrow_reads_discarded > 0:
    logging.info('\t%d out of %d escrow reads had to be discarded because their singular weight was too low, but %d (%7.3f%%) were still retained',
                Reads.total_escrow_reads_discarded,
                Reads.total_escrow_reads,
                Reads.total_escrow_reads_kept,
                Reads.total_escrow_reads_kept * 100 / Reads.reads_w_perfect_alignments)
    
if singular_fraud_alert or escrow_fraud_alert:
    logging.critical('[FATAL] There were signs of voter fraud. See the votes file (%s) for the complete picture.', vote_file)
    logging.critical("\tVoter fraud committed by " +str(num_singular_fraud_reads)+ " out of " +str(Reads.total_singular_reads_after_recruitment)+ " singular reads in pairs")
    logging.critical("\tVoter fraud committed by " +str(num_fraud_reads)+ " out of " +str(Reads.total_escrow_reads)+ " escrow reads")
    sys.exit(1)

###############################################################################
# Sanity Check
###############################################################################
# for name, strain in all_strain_obj.items():
#     singular_vote_ratio = escrow_vote_ratio = 0
#     if strain.num_singular_reads > 0:
#         singular_vote_ratio = strain.cum_primary_votes/strain.num_singular_reads
#     if strain.num_escrow_reads > 0:
#         escrow_vote_ratio = strain.cum_escrow_votes/strain.num_escrow_reads
                             
#     cumulative_vote = (strain.cum_escrow_votes + strain.cum_primary_votes)/Reads.total_reads_aligned
#     if cumulative_vote > 10e-2:
# #     if strain.num_singular_reads > 0:
#         logfile.info(name + '\t' + str(strain.num_singular_reads) + '\t' + str(strain.num_escrow_reads) + '\t' +  str(singular_vote_ratio) + '\t' + str(escrow_vote_ratio) + '\t' + str(cumulative_vote))

###############################################################################
# Calculate unique number of bases covered for each genome
# Create strain stats file
# Create relative abundance output file
###############################################################################

logging.info('Computing depth and coverage for each strain in the database based on escrow alignments ...')

abundance_df = pd.DataFrame()
stats_df = pd.DataFrame()

# pysam.sort('-o','tmp.bam', escrowBamfile)
# os.rename('tmp.bam',escrowBamfile)
# pysam.index(escrowBamfile)
sorted_escrowBamfile = sort_and_index(escrowBamfile)

escrow_bed_coverage = Strains.calculate_coverage(sorted_escrowBamfile)
i = 0
escrow_depth_file = prefix +'.ninjaMap.escrow_depth.csv'
with open(escrow_depth_file, 'w+') as ed_handle:
    ed_handle.write(f"strain\tcoverage\tdepth\tdepth_sd\tdepth_coeff_variance\n")
    for name, strain in all_strain_obj.items():
        i += 1
        logging.info("\t[%d/%d]Searching for escrow support for :\t%s",i, Strains.total_strains, name)
        (coverage, wt_depth, depth_sd, depth_variance) = strain.calculate_genome_coverage(escrow_bed_coverage, singular = False)
        if coverage > 0:
            logging.info(f"\t... Calculated Escrow Coverage = {coverage}%; Weighted Depth = {wt_depth}x; Depth StDev = {depth_sd}; Variation in Depth = {depth_variance}%")
            ed_handle.write(f"{strain.name}\t{coverage}\t{wt_depth}\t{depth_sd}\t{depth_variance}\n")
        my_frac = strain.calculate_read_fraction()
        logging.info(f"\t... Calculated relative abundance : {my_frac}%")
        strain_stats_df = strain.compile_general_stats()
        if strain_stats_df is not None:
            stats_df = pd.DataFrame.add(stats_df, strain_stats_df, fill_value = 0)
        strain_abundance_df = strain.compile_by_abundance()
        if strain_abundance_df is not None:
            abundance_df = pd.DataFrame.add(abundance_df, strain_abundance_df, fill_value = 0)

logging.info('Writing strain stats file ...')
stats_df.to_csv(strain_stats_file, index_label='Strain_Name')

logging.info('Creating the relative abundance file ...')
abundance_df.to_csv(abundance_output_file, index_label='Strain_Name')

###############################################################################
# Stats
###############################################################################
stats = open(stats_file, 'w')
# Header
stats.write('File_Name,Reads_Aligned,Reads_wPerfect_Aln,Reads_wSingular_Votes,Reads_wEscrowed_Votes,Discarded_Reads_w_Perfect_Aln\n')
stats.write(
    prefix +','+
    str(Reads.total_reads_aligned) +','+ 
    str(Reads.reads_w_perfect_alignments) +','+
    str(Reads.total_singular_reads_after_recruitment) +','+
    str(Reads.total_escrow_reads_kept) +','+
    str(Reads.total_escrow_reads_discarded) +'\n')
stats.close()

###############################################################################
# Running an independent voting check
###############################################################################
logging.info('Running an independent voter fraud check ...')
fraud_committed = Reads.check_for_fraud(vote_file, Reads.total_escrow_reads_discarded)
if fraud_committed:
    logging.critical('\tEvidence of voter fraud was found! :(')
else:
    logging.info('\tNo evidence of voter fraud detected! :)')

###############################################################################
# The End
###############################################################################
end = timer()
logging.info('Completed in %s (d:h:m:s)', human_time(end - start))
