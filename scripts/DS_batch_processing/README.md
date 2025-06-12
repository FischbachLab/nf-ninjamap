# How to bulk generate seedfiles and parameter files for all of the DS samples?

### &bull; The metadata format: a 6-column tsv or cvs file with the following 6 headers (db_file_path is optional)
`SampleID	R1	R2	DB_name Project	Output_path`

[View example metadata file content](https://github.com/FischbachLab/nf-ninjamap/blob/main/scripts/DS_batch_processing/examples/DS_metadata.tsv)

### 0. Install python library (optional)
```{bash}
pip install pandas boto3
```
### 1. Generate the seedfiles, parameters and batch files for all DS samples, automatically 
```{bash}
python create_DS_ninjamap_jobs.py --input DS_metadata.tsv --output_dir ./outfiles_dir --batch_name my_batch_file_name.sh
```
### 2. Run the following command to submit batch jobs, once the script finishes successfully
```{bash}
bash  ./outfiles_dir/my_batch_file_name.sh
```
#### &bull; An example of the generated batch job in ./outfiles_dir/my_batch_file_name.sh
```{bash}
aws batch submit-job \
          --job-name nf-ninjamap-20250418_MITI-001_DS072-mNGS_HVFW7BGYW_Q30-L50 \
          --job-queue priority-maf-pipelines \
          --job-definition nextflow-production \
          --container-overrides command="fischbachlab/nf-ninjamap, \
         "-params-file", "s3://genomics-workflow-core/Results/Ninjamap/MITI-001v3_20240604/bespoke-DS/parameter_files/20250418_MITI-001_DS072-mNGS_HVFW7BGYW_Q30-L50.parameters.json" " 
```
