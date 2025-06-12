#!/usr/bin/env python3
"""Automatically generate seedfiles, parameters and batch files for DS samples for ninjaMap pipeline,
   then upload to s3 bucket
USAGE: python create_DS_ninjamap_jobs.py --input DS_metadata.tsv --output_dir ./outfiles_dir --batch_name my_batch_file_name.sh

Once the script is done successfully, run the following command
###################################################
bash  ./outfiles_dir/my_batch_file_name.sh
##################################################

# Install python library
# pip install pandas boto3
"""

import boto3
from botocore.exceptions import ClientError
import logging
import os, sys
import pandas as pd
#import numpy as np
import logging
import argparse
#from io import StringIO
import json

MITI_NINJAMAP_INDEX="s3://maf-versioned/ninjamap/Index/MITI-001-DS"
#s3_output_path="s3://genomics-workflow-core/Results/Ninjamap/test"

def usage():
    usage = """
    python create_DS_ninjamap_jobs.py \\
        --input input_metadata_sheet.tsv  input metadata sheet in tsv or csv format with headers \\
        --output_dir  /local/path/to/output/dir  path to the output batch scripts \\
        --batch_name batch_file_name.sh default: batch_MITI_DS_ninjamap.sh
    """

    # Making default argument list structures
    p = argparse.ArgumentParser(usage=usage)
    # Required
    p.add_argument(
        "-i",
        "--input",
        dest="input",
        action="store",
        type=str,
        required=True,
        help="The path to input metadata sheet in tsv or csv format with headers",
    )
    p.add_argument(
        "-o",
        "--output_dir",
        dest="output",
        action="store",
        type=str,
        required=True,
        help="The local output directory name",
    )
    p.add_argument(
        "-n",
        "--batch_name",
        dest="batch_name",
        action="store",
        type=str,
        default="batch_MITI_DS_ninjamap.sh",
        required=False,
        help="The output aws batch file name",
    )
  
    return vars(p.parse_args())
"""
  
"""
def upload_to_s3(file, s3_path, bucket_name):
    """Upload a file to an S3 bucket
    :param file_name: File to upload
    :param bucket: Bucket to upload to
    :param s3_path: S3 path
    :param project: project name
    :return: s3 file path
    """
    #if object_name is None:
    #    object_name = os.path.basename(file_name)

    client = boto3.client('s3')
    try:
        client.upload_file(file, bucket_name, s3_path)
        #client.put_object(Bucket=bucket_name, Key=s3_path, Body=buffer.getvalue())
    except ClientError as e:
        logging.error(e)
        #return False
    print(f"File upload successfully to s3://{bucket_name}/{s3_path}")
    return f"s3://{bucket_name}/{s3_path}"


def download_file_from_s3( bucket_name, key_prefix, destination_prefix, aws_profile=None):
    """
    transfer a file to/from AWS S3
    """
    if destination_prefix == ".":
        destination_prefix = os.path.curdir
    destination = os.path.join(destination_prefix, os.path.basename(key_prefix))
    # logging.info(f"Downloading {key_prefix}")
    if aws_profile is not None:
        # logging.info(f"Downloading using profile {aws_profile}")
        s3 = boto3.session.Session(profile_name=aws_profile).client("s3")
    else:
        # logging.info(f"Downloading using default profile")
        s3 = boto3.client("s3")
    try:
        s3.download_file(bucket_name, key_prefix, destination)
    except botocore.exceptions.ClientError as b:
        logging.error(
            f"Failed to retrieve object with bucket_name:{bucket_name}; key_prefix:{key_prefix}; destination:{destination}"
        )
    return destination

def verify_s3_path(s3path):
    """
    take s3 path
    return nothing
    """
    path_list = s3path.replace("s3://", "").split("/")
    bucket = path_list.pop(0)
    s3_key = "/".join(path_list)
    #return (bucket, obj_key)

    client = boto3.client('s3')
    content = client.head_object(Bucket=bucket,Key=s3_key)
    if content.get('ResponseMetadata',None) is not None:
        print (f"File exists - s3://{bucket}/{s3_key} ")
    else:
        print (f"File does not exist - s3://{bucket}/{s3_key} ")
    #return 

def create_seedfile(sample_name, R1, R2, DB, project, s3_output_path, batch_file_name, out_dir):
    print (f"{sample_name}")
    print (f"{R1}")
    print (f"{DB}")
    # test s3 path 
    verify_s3_path(R1)
    verify_s3_path(R2)
    sd = pd.DataFrame(
        [
            {
                "sampleName": str(sample_name),
                "R1": str(R1),
                "R2": str(R2),
            }   
        ]
    )
    
    s3_output_path= s3_output_path.rstrip("/")  # Removes trailing slash if present
    print(s3_output_path)  

    seedfile = f"{project}.seedfile.csv"
    path_list = s3_output_path.replace("s3://", "").split("/")
    bucket = path_list.pop(0)
    s3_key = "/".join(path_list)

      # Convert DataFrame to CSV
    #csv_buffer = StringIO()
    #sd.to_csv(csv_buffer, index=False)
    csv_str= sd.to_csv(index=False) 

    print(f"{bucket}")
    print (f"{s3_key}/{seedfile}")

      # Write to a file
    with open(f"{out_dir}/DS_seedfiles/{seedfile}", "w", encoding="utf-8") as file:
        file.write(csv_str)

    seedfile_s3_path = upload_to_s3(f"{out_dir}/DS_seedfiles/{seedfile}", f"{s3_key}/seedfiles/{seedfile}", bucket)
    #print (f"{seedfile_s3_path}")
    parameter_s3_path = create_parameter_file(s3_output_path, seedfile_s3_path, DB, project, out_dir)
    #print (f"{parameter_s3_path}")
    create_batch_file(parameter_s3_path, project, out_dir, batch_file_name)
    

def create_parameter_file(s3_output_path, seedfile_file_path, db_name, project, out_dir):

    parameter_df = pd.DataFrame(
            {
                "seedfile": seedfile_file_path,
                "project": project,
                "db_path": MITI_NINJAMAP_INDEX,
                "db": db_name,
                "db_prefix": db_name,
                "minQuality": 30,
                "minLength": 50,
                "coreNum": 16,
                "sampleRate": 1,
                "coverage": 1,
                "debug": 1,
                "output_path": s3_output_path
            }, index=[0]  
    )

    #json_buffer = StringIO()
    #parameter_df.to_json(json_buffer, orient="records", index=False)  # "records" format makes it a list of dictionaries

    json_list = parameter_df.to_json(orient="records") 
    # Remove the outer list brackets and write to a file
    json_objects = json.loads(json_list)  # Convert JSON string to a list of dictionaries
    # Pretty-print JSON and ensure S3 paths are displayed correctly
    json_str = "\n".join(json.dumps(obj, indent=4, ensure_ascii=False) for obj in json_objects)  # Format without brackets

    #json_pretty = json.dumps(json.loads(json_file), indent=4, ensure_ascii=False)

    # parameter file name
    parameterfile = f"{project}.parameters.json"
    s3_path_list = s3_output_path.replace("s3://", "").split("/")
    bucket = s3_path_list.pop(0)
    s3_key = "/".join(s3_path_list)

     # Write to a file
    with open(f"{out_dir}/DS_parameter_files/{parameterfile}", "w", encoding="utf-8") as file:
        file.write(json_str)

    parameter_s3_path = upload_to_s3(f"{out_dir}/DS_parameter_files/{parameterfile}", f"{s3_key}/parameter_files/{parameterfile}", bucket)

    return parameter_s3_path 


def create_batch_file(paramter_file_path, project, out_dir, batch_file_name):
    """
    print the following batch command lines into a file and wait for 5s
    aws batch submit-job \
        --job-name nf-ninjamap-MITI \
        --job-queue priority-maf-pipelines \
        --job-definition nextflow-production \
        --container-overrides command="fischbachlab/nf-ninjamap, \
        "-params-file", "s3://genomics-workflow-core/Results/Ninjamap/parameters/example_parameters.json" " 
    """
    cmd_file=batch_file_name
    batch_cmd = f"""aws batch submit-job \\
          --job-name nf-ninjamap-{project} \\
          --job-queue priority-maf-pipelines \\
          --job-definition nextflow-production \\
          --container-overrides command=\"fischbachlab/nf-ninjamap, \\
         \"-params-file\", \"{paramter_file_path}\" " 

sleep 5

"""
    #batch_cmd.to_csv(out_dir/cmd_file, header = False, index = False)
    with open( f"{out_dir}/{cmd_file}", "a", encoding="utf-8") as file:
        file.write(batch_cmd)

    print(f"Batch script {out_dir}/{cmd_file} written successfully!")
    

def process_input_metadata2(infile):
    """Automatically detects if a file is CSV or TSV and loads it into a DataFrame."""
    with open(infile, 'r', encoding='utf-8') as f:
        sample = f.readline()  # Read the first line
    
    if '\t' in sample:
        delimiter = '\t'  # Likely a TSV file
    elif ',' in sample:
        delimiter = ','  # Likely a CSV file
    else:
        raise ValueError("Unable to detect file format. Ensure metadata is a CSV or TSV.")

    df = pd.read_csv(infile, delimiter=delimiter)
    return df

# Not working 
def process_input_metadata(infile):
    try:
    #read tsv file
        df = pd.read_table(infile) #, names=['sampleID', 'R1', 'R2', 'DB_name', 'Project', 'Output_path'])
    except (pd.errors.ParserError, FileNotFoundError):
        df = pd.read_csv(infile) #, names=['sampleID', 'R1', 'R2', 'DB_name', 'Project', 'Output_path'])
    return df
    
def remove_file_if_exists(file_path):
    """Check if a file exists and remove it."""
    if os.path.exists(file_path):  # Check if the file exists
        os.remove(file_path)  # Remove the file
        print(f"Old file '{file_path}' removed successfully!")
    #else:
    #    print(f"File '{file_path}' does not exist.")


def main():
    args = usage()
    input_file = args["input"]
    batch_name = args["batch_name"]
    script_outdir = args["output"]

    #s3_output_path="s3://genomics-workflow-core/Results/Ninjamap/test"
    os.makedirs(script_outdir , exist_ok=True)
    os.makedirs(f"{script_outdir}/DS_parameter_files", exist_ok=True)
    os.makedirs(f"{script_outdir}/DS_seedfiles", exist_ok=True)

    remove_file_if_exists(f"{script_outdir}/{batch_name}")

    db_df = pd.DataFrame()
    # Read input
    db_df = process_input_metadata2(input_file)

    print(db_df.head(5))
    db_df.apply(lambda row: create_seedfile(sample_name = row['SampleID'], 
                                            R1 = row['R1'], 
                                            R2 = row['R2'],
                                            DB = row['DB_name'],
                                            project = row['Project'],
                                            s3_output_path = row['Output_path'], 
                                            batch_file_name = batch_name,
                                            out_dir = script_outdir),                                      
                                            axis=1).dropna(axis='index')

if __name__ == '__main__':
     main()

# Usage
#upload_file('path/to/file.txt', 'mybucket', 'file.txt')