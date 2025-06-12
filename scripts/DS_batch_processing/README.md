# How to bulk generate seedfiles and parameter files for all of the DS samples?

### &bull; The metadata format: a four-colum tsv or cvs file with the following 6 headers (db_file_path is optional)
`SampleID	R1	R2	DB_name Project	Output_path`
[View example metadata file content](https://github.com/FischbachLab/nf-ninjamap/blob/main/scripts/DS_batch_processing/examples/DS_metadata.tsv)

### &bull; Install python library (optional)
```{bash}
pip install pandas boto3
```
### &bull; Generate the seedfile, parameter and batch files for DS samples for ninjaMap pipeline, automatically 
```{bash}
python create_DS_ninjamap_jobs.py --input DS_metadata.tsv --output_dir ./outfiles_dir --batch_name my_batch_file_name.sh
```
### &bull; Once the script is generated successfully, run the following command
```{bash}
bash  ./outfiles_dir/my_batch_file_name.sh
```
#### &bull; An example of generated batch job by the above script 
```{bash}
aws batch submit-job \
          --job-name nf-ninjamap-20250418_MITI-001_DS072-mNGS_HVFW7BGYW_Q30-L50 \
          --job-queue priority-maf-pipelines \
          --job-definition nextflow-production \
          --container-overrides command="fischbachlab/nf-ninjamap, \
         "-params-file", "s3://genomics-workflow-core/Results/Ninjamap/MITI-001v3_20240604/bespoke-DS/parameter_files/20250418_MITI-001_DS072-mNGS_HVFW7BGYW_Q30-L50.parameters.json" " 
```
