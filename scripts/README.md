Aggregated output files for each study
====================


### The following steps show how to aggregate the ninjamap results from individual samples stored on s3 into the summary tables and plots. 
#### The s3 path should be in the format `<s3_nm_base_path>/<db_name>/<study_name>`

1. Update the paths and variables in the following sctipt, then run the R Markdown file in RStudio
```{bash}
scripts/aggregation/aggregate_results_template.Rmd
```
2. Update the nm_basedir path, which is the library path to aggregate_data_util_functions.R
3. Update s3_nm_base_path, e.g., s3://genomics-workflow-core/Results/Ninjamap  
4. Update db_name. e.g., hCom2
5. Update study_name, e.g., hCom2-DietaryIntervention
6. Update the output dir: workdir

Outputs files
====================

### The aggregated output files are organized into 6 files.

1. **\*.covDepth.csv**: this file shows the average coverage depth per strain by samples.
2. **\*.host_contaminants.csv**: this file shows the detected host contaminants (Human and Mouse) by samples if the unaligned read rate is over 5%.
3. **\*.long.csv**: this is the long format of three files (\*.readFraction.csv, \*.covDepth.csv and \*.percCoverage.csv)
4. **\*.percCoverage.csv**: this file shows the average coverage per strain in percentage by samples.
5. **\*.reads_stats.csv**: this file shows the reads statistics in read numbers by samples. (Reads_wPerfect_Aln (Col.G) / Reads_Aligned(Col.F) is equal to the sum of all abundance numbers for a sample in readFraction.)
6. **\*.readFraction.csv**: this file shows the relative abundance in the defined community in percentage by samples. (The relative abundance is denoted by the fraction of the total QCed filtered reads in a sample that aligned perfectly to at least one strain in the database. Since that almost never happens, the fractions don’t sum to 100%.)

### The ninjaMap read fate definitions in waffle plot.

1. **Primary**: reads that align to a single reference within the queried database
2. **Escrow**: reads that align to multiple references within the queried database. Escrow reads are divided based on the Primary read fraction of the relevant references.
3. **QC Fail**: reads whose length < 50 after read trimming at a threshold of Q<30 by default
4. **Missed**: reads that have some amount of mismatch against the queried db but still align to a targeted reference within the queried db ( up to 18% mismatch)
5. **Unaligned**: reads that do not align to the queried db (designation from bowtie2)
6. **Discarded**: reads that perfectly match too many strains in the queried db such that they cannot be distributed among the db strains. These reads represent coverage of highly conserved regions.