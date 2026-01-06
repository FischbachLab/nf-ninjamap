#!/bin/bash -x
# shellcheck disable=SC2086
# shellcheck disable=SC2154

#set -e
set -u
set -o pipefail

in_abundance=$1
singular_bam=$2
escrow_bam=$3

bedN="${bedN:-100000}" # bedfile threshold 0.1M 

NINJA_OUTPUT="."
GENOME_COV_OUTPUT="."

if [ ${coverage} -eq 1 ]; then
    echo "Compute Singular coverage and depth"
    coreN=$(grep -c ^processor /proc/cpuinfo)
    samtools faidx ${REF} && cut -f1,2 ${REF}.fai > Ninja.genomes
    cut -f1  Ninja.genomes | sed  's/_Node.*//' | sort | uniq  > DBGenomeNameList.txt


    s="singular"
    samtools sort -@ ${coreN} ${singular_bam} > ${SAMPLE_NAME}.singular_sorted.bam
    bedtools bamtobed -i ${SAMPLE_NAME}.singular_sorted.bam | cut -f1,2,3 > ${SAMPLE_NAME}.${s}.bed
    samtools index  ${SAMPLE_NAME}.singular_sorted.bam
    printf "Strain_Name\tSingular_Coverage\tSingular_Depth\tSingular_Bases\n" > ${s}_summary_depth.tsv
    
    for i in $(cat DBGenomeNameList.txt)
    do
        if [ $(grep -c "^$i" ${SAMPLE_NAME}.${s}.bed) -eq 0 ];
        then
          > tmp.bed
        else
          grep "^$i" ${SAMPLE_NAME}.${s}.bed | sort | uniq | cut -f1,2,3  > tmp.bed
        fi
        grep "^$i"  Ninja.genomes > tmp.genome
        bedtools genomecov -i tmp.bed -g tmp.genome -bga > tmps.bedg
        
        zero=$(awk '$4==0 {bpCountZero+=($3-$2)} END {print bpCountZero}' tmps.bedg)
        nonzero=$(awk '$4>0 {bpCountNonZero+=($3-$2)} END {print bpCountNonZero}' tmps.bedg)
        if [ -z $zero ];  then zero=0; fi
        if [ -z $nonzero ];  then nonzero=0; fi
        cov_pct=$(bc <<< "scale=4; ($nonzero / ($zero + $nonzero))*100")
        
        totalbases=$(awk '{t+=$4*($3-$2)} END {print t}' tmps.bedg)
        if [ $nonzero -eq 0 ] || [ $totalbases -eq 0 ] ;
        then
           cov_depth=0
        else
           cov_depth=$(bc <<< "scale=4; ($totalbases / $nonzero)")
        fi
        printf "${i}\t${cov_pct}\t${cov_depth}\t${nonzero}\n" >> ${s}_summary_depth.tsv
     done
     #cut -f2,3,4 ${s}_summary_depth.tsv | paste ${in_abundance} -  | awk '{print $1","$2","$3","$4}' > tmp.ninjaMap.abundance.csv
    merge_abundance_coverage_table.py ${SAMPLE_NAME}  ${in_abundance}  ${s}_summary_depth.tsv

     echo "Compute Escrow coverage and depth"
     s="escrow"
     samtools sort -@ ${coreN} ${escrow_bam} > ${SAMPLE_NAME}.escrow_sorted.bam
     bedtools bamtobed -i ${NINJA_OUTPUT}/${SAMPLE_NAME}.escrow_sorted.bam  | cut -f1,2,3 > ${NINJA_OUTPUT}/${SAMPLE_NAME}.${s}.bed
     samtools index ${NINJA_OUTPUT}/${SAMPLE_NAME}.escrow_sorted.bam 
     printf "Strain_Name\tEscrow_Coverage\tEscrow_Depth\tEscorw_Bases\n" > ${NINJA_OUTPUT}/${s}_summary_depth.tsv
     for i in $(cat ${GENOME_COV_OUTPUT}/DBGenomeNameList.txt)
     do
         if [ $(grep -c "^$i" ${NINJA_OUTPUT}/${SAMPLE_NAME}.${s}.bed) -eq 0 ];
         then
            > ${GENOME_COV_OUTPUT}/tmp.bed
         else
           grep "^$i" ${NINJA_OUTPUT}/${SAMPLE_NAME}.${s}.bed | sort | uniq | cut -f1,2,3  > ${GENOME_COV_OUTPUT}/tmp.bed
         fi
         grep "^$i"  ${GENOME_COV_OUTPUT}/Ninja.genomes > ${GENOME_COV_OUTPUT}/tmp.genome
         bedtools genomecov -i ${GENOME_COV_OUTPUT}/tmp.bed -g ${GENOME_COV_OUTPUT}/tmp.genome -bga > ${GENOME_COV_OUTPUT}/tmps.bedg
         #compute coverage
         zero=$(awk '$4==0 {bpCountZero+=($3-$2)} END {print bpCountZero}' ${GENOME_COV_OUTPUT}/tmps.bedg)
         nonzero=$(awk '$4>0 {bpCountNonZero+=($3-$2)} END {print bpCountNonZero}' ${GENOME_COV_OUTPUT}/tmps.bedg)
         if [ -z $zero ];  then zero=0; fi
         if [ -z $nonzero ];  then nonzero=0; fi

         cov_pct=$(bc <<< "scale=4; ($nonzero / ($zero + $nonzero))*100")
        
         # compute depth
         totalbases=$(awk '{t+=$4*($3-$2)} END {print t}' ${GENOME_COV_OUTPUT}/tmps.bedg)
         if  [ $nonzero -eq 0 ] || [ $totalbases -eq 0 ];
         then
            cov_depth=0
         else
            cov_depth=$(bc <<< "scale=4; ($totalbases / $nonzero)")
         fi
         printf "${i}\t${cov_pct}\t${cov_depth}\t${nonzero}\n" >> ${NINJA_OUTPUT}/${s}_summary_depth.tsv
      done

      merge_abundance_coverage_table.py ${SAMPLE_NAME} ${SAMPLE_NAME}.ninjaMap.abundance.csv ${NINJA_OUTPUT}/${s}_summary_depth.tsv
      #cut -f2,3,4 ${NINJA_OUTPUT}/${s}_summary_depth.tsv | paste ${NINJA_OUTPUT}/tmp.ninjaMap.abundance.csv -  | awk '{print $1","$2","$3","$4}' > ${NINJA_OUTPUT}/tmp2.ninjaMap.abundance.csv
      #mv ${NINJA_OUTPUT}/tmp2.ninjaMap.abundance.csv ${NINJA_OUTPUT}/${SAMPLE_NAME}.ninjaMap.abundance.csv


    if [ ${debug} -eq 1 ];then
        echo "Enter into the debug mode"
        mkdir -p "${NINJA_OUTPUT}/debug/all"
        mv ${NINJA_OUTPUT}/*_sorted.bam "${NINJA_OUTPUT}/debug/all"
        mv ${NINJA_OUTPUT}/*_sorted.bam.bai "${NINJA_OUTPUT}/debug/all"
        mv ${NINJA_OUTPUT}/*.bed "${NINJA_OUTPUT}/debug/all"

        for s in {singular,escrow}; 
        do
            mkdir -p "${NINJA_OUTPUT}/debug/${s}"
            bed_line=0
            for i in $(cat ${GENOME_COV_OUTPUT}/DBGenomeNameList.txt)
            do
                bed_line=`grep -c "^$i" ${NINJA_OUTPUT}/debug/all/${SAMPLE_NAME}.${s}.bed`
                if [ ${bed_line} -eq 0 ];
                then
                    > ${GENOME_COV_OUTPUT}/tmp.bed
                    else
                    grep "^$i" ${NINJA_OUTPUT}/debug/all/${SAMPLE_NAME}.${s}.bed | sort | uniq | cut -f1,2,3  > "${NINJA_OUTPUT}/debug/${s}/${i}.bed"
                    
                    if [ ${bed_line} -lt ${bedN} ];
                    then
                        
                        bedtools intersect -abam ${NINJA_OUTPUT}/debug/all/${SAMPLE_NAME}.${s}_sorted.bam -b "${NINJA_OUTPUT}/debug/${s}/${i}.bed" > "${NINJA_OUTPUT}/debug/${s}/${i}_${s}.bam"
                        samtools index -@ ${coreN} "${NINJA_OUTPUT}/debug/${s}/${i}_${s}.bam"
                    fi
                fi
            done
        done    
    fi
fi

rm ${GENOME_COV_OUTPUT}/tmp.bed