#############################################################################
# 
# Run sanity checks on fastq files, logging md5sum, sha256, and file_size of
# local fastq files and on S3. There is no retrieval of metadata from an url.
#
# Sample commandline:
# checkfiles_lite.py --out outfile.txt --err err_log.txt --s3-bucket bucket_name 
#       --s3-dir s3_dir  --local-dir ec2_local_dir --hca-manifest hca-manifest.tsv
#
##############################################################################

# Install slackclient==1.3.1
import datetime
import time
import os.path
import sys
import re
import boto3
import pandas as pd
from os import listdir
import subprocess

# Entries to be made in the result dictionary
features = ['local_md5','local_sha','hca_size','hca_sha256','s3_file_size']

# Calculate md5sum and sha256 using subprocess
def calc_hash(result, local_file, err):
    md5_encoded = subprocess.check_output(['md5sum {}'.format(local_file)],
            shell=True,
            executable='/bin/bash',
            stderr=subprocess.STDOUT)
    md5 = md5_encoded.decode('utf-8').split(' ')[0]
    sha_encoded = subprocess.check_output(['sha256sum {}'.format(local_file)],
            shell=True,
            executable='/bin/bash',
            stderr=subprocess.STDOUT)
    sha = sha_encoded.decode('utf-8').split('  ')[0]
    
    # Add new entry to result if needed, add md5 and sha info
    local_file_split = local_file.split("/")
    local_filename = local_file_split[len(local_file_split)-1]
    if local_filename not in result:
        result[local_filename] = {}
        err.write("local file unique: {}".format(local_file))
        err.flush()
    result[local_filename]['local_md5'] = md5
    result[local_filename]['local_sha'] = sha


# Retrieve file_sha256 and file_size
def check_manifest(result, hca_manifest, err):
    manifest = pd.read_csv(hca_manifest, sep="\t", header=0)
    for index, row in manifest.iterrows():
        if row['file_name'] not in result:
            result[row['file_name']] = {}
            err.write("manifest file unique: {}".format(row['file_name']))
            err.flush()
        result[row['file_name']]['hca_size'] = row['file_size']
        result[row['file_name']]['hca_sha256'] = row['file_sha256']
    
    return result


# Retreive e_tag and file size for s3_object and store in result
def check_s3_file(result, s3_object, s3_dir=None):
   
    # If there is a subdirectory in s3, remove from key to obtain filename
    if s3_dir:
        key = s3_object.key.split("/")
        filename = key[len(key)-1]
    else:
        filename = str(s3_object.key)
    
    # Add new entry to result if needed, add file_size and etag info 
    if filename not in result:
        result[filename] ={}
    result[filename]['s3_file_size'] = s3_object.size
    
    return result 


# Go through s3, local_dir, and hca_manifest to obtain and/or calculate file metadata
def run(out, err, s3_bucket, local_dir, hca_manifest, s3_dir=None):
    result = {}

    # access s3 bucket, and iterate through files in bucket
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(s3_bucket)
    if s3_dir:
        bucket_files = list(bucket.objects.filter(Prefix=s3_dir))
    else:
        bucket_files = bucket.objects.all()

    # For each object in bucket, if it is a fastq file, fetch information and run check_s3_file
    for s3_object in bucket_files:
        if s3_object.key.endswith('fastq.gz'):
            check_s3_file(result, s3_object, s3_dir)

    # Go through hca manifest to retrieve sha25 and file_size
    check_manifest(result, hca_manifest, err)

    # calculate md5sum and sha256 on local files
    local_files = os.listdir(local_dir)
    for local_file in local_files:
        calc_hash(result, local_dir+local_file, err)

    # Make sure all features are in result, if not, default to ''
    for f in result:
        for feat in features:
            if feat not in result[f]:
                result[f][feat] = ''

    # Print to log
    out.write("file_name\ts3_file_size\thca_size\thca_sha256\tlocal_sha\tlocal_md5\n")
    for filename in result:
        out.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(filename, result[filename]['s3_file_size'], result[filename]['hca_size'],
                result[filename]['hca_sha256'], result[filename]['local_sha'], result[filename]['local_md5']))

    out.flush()
    output_filename = out.name
    out.close()
    error_filename = err.name
    err.close()


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Update file status", 
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--out', '-o', type=argparse.FileType('w'), default=sys.stdout,
        help="file to write json lines of results with or without errors")
    parser.add_argument(
        '--err', '-e', type=argparse.FileType('w'), default=sys.stderr,
        help="file to write json lines of results with errors")
    parser.add_argument(
        '--s3-bucket', required=True,
        help="s3 bucket name that contains the fastq files")
    parser.add_argument(
        '--local-dir', required=True,
        help="local directory that contains the fastq files")
    parser.add_argument(
        '--hca-manifest', type=argparse.FileType('r'), required=True,
        help="hca manifest tsv file from the DCP for fastq files")
    parser.add_argument(
        '--s3-dir', default='',
        help="directory on the S3 where fastq files are located")

    args = parser.parse_args()
    run(**vars(args))


if __name__ == '__main__':
    main()
