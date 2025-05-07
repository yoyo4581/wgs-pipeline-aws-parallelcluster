import json
import urllib.parse
import boto3
import os
import time
import csv
import io
import yaml

print('Submitting to EC2 Login...')

ssm_client = boto3.client('ssm')
ec2_client = boto3.client('ec2')

dynamodb = boto3.resource('dynamodb')
table = dynamodb.Table('SRATable')

ttl_seconds = int(time.time()) + 7*24*60*60


snakemake_source = 's3://wgs-snakemake-files-yoyo458/'
efs_mount_path = '/home/ec2-user/snakemake_files/'

# finds the instance id of the head node
def get_instance_id_by_name(instance_name):
    response = ec2_client.describe_instances(
        Filters=[
            {
                'Name': 'tag:Name',
                'Values': [instance_name]
            }
        ]
    )
    
    if len(response['Reservations']) > 0:
        instances = response['Reservations'][0]['Instances']
        if len(instances) > 0:
            return instances[0]['InstanceId']
    else:
        raise Exception(f"No instances found with name: {instance_name}")


# sync snakemake workflow files from s3 to efs mounted directory on the head node
def update_snakefiles(instance_id, s3_path):
    document_name = 'AWS-RunShellScript'
    efs_mount_path = '/workflows/snakemake_files'

    commands = [
        f'mkdir -p {efs_mount_path}',
        f'aws s3 sync {s3_path} {efs_mount_path}',
        f'sudo su - ec2-user',
    ]

    response = ssm_client.send_command(
        InstanceIds=[instance_id],
        DocumentName=document_name,
        Parameters={'commands': commands},
        CloudWatchOutputConfig={
            'CloudWatchLogGroupName': '/aws/ssm/SnakefileSync',
            'CloudWatchOutputEnabled': True
        }
    )

    command_id = response['Command']['CommandId']
    print(f"Command ID: {command_id}")
    return command_id


def initiate_snakemake_slurm():
    # Define tag to search for
    instance_name = 'HeadNode'

    try:
        # Get the instance ID by name
        instance_id = get_instance_id_by_name(instance_name)

        print(f"Instance ID: {instance_id}")
        command_response = update_snakefiles(instance_id, snakemake_source)

    except Exception as e:
        print(f"Error: {str(e)}")
        raise
    
    return instance_id


