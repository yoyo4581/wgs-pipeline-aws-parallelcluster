import json
import urllib.parse
import boto3
import os
import time
import csv
import io
import yaml
from task_submit import initiate_snakemake_slurm

print('Loading function')

# Initialize the AWS services
ecs_client = boto3.client('ecs')
s3_client = boto3.client('s3')
secret_client = boto3.client('ssm')
ec2_client = boto3.client('ec2')
dynamodb = boto3.resource('dynamodb')
table = dynamodb.Table('SRATable')



ttl_seconds = int(time.time()) + 7*24*60*60


def generate_config(queued_entries):

    config = {'metadata': {}}

    for entry in queued_entries:

        experiment_id = entry['ExperimentId']
        srr_id = entry['SRRId']

        if experiment_id in config['metadata']:
            config['metadata'][experiment_id].append(srr_id)
        else:
            config['metadata'][experiment_id] = [srr_id]
    
    return config

# Acquires an s3 bucket event finds the bucket and key, and processes csv file and stores sample entry in a dynamodb table.
# Also generates a config.yaml file and stores it in the s3 bucket wgs-genomics-files/input this will be used to structure workflow subdirectorie in snakemake
# After storing initiates a snakemake job on the HPC cluster using slurm
def lambda_handler(event, context):
    print("Received event: " + json.dumps(event, indent=2))

    # Get the bucket and object key from the S3 event
    bucket = event['Records'][0]['s3']['bucket']['name']
    key = urllib.parse.unquote_plus(event['Records'][0]['s3']['object']['key'], encoding='utf-8')
    
    try:
        # Get the Snakefile from S3
        response = s3_client.get_object(Bucket=bucket, Key=key)
        print("CONTENT TYPE: " + response['ContentType'])

        csv_content = list(csv.reader(io.StringIO(response['Body'].read().decode('utf-8'))))
        header = csv_content[0]
        
        fields = {header[i]: i for i in range(len(header))}
        print(len(fields))

        for entry in csv_content[1:]:
            print(len(entry))
            item = {
                'SRRId': entry[fields['Run']],
                'ExperimentId': entry[fields['Experiment']],
                'CellLine': entry[fields['cell_line']],
                'MetadataBlob': json.dumps({
                    'LibraryLayout': entry[fields['LibraryLayout']],
                    'AvgSpotLength': entry[fields['AvgSpotLen']],
                    'Platform': entry[fields['Platform']],
                }),
                'ProcessingStatus': 'PENDING',
                'TimeToExist': ttl_seconds
            }

            try:
                # Put the item into the table
                response = table.put_item(Item=item, ConditionExpression='attribute_not_exists(SRRId)')
                print("DynamoDb Put Item Response: ", json.dumps(response, indent=2))
            except Exception as e:
                if e.response['Error']['Code'] == 'ConditionalCheckFailedException':
                    print(f"Item with SRRId {item['SRRId']} already exists!")
                else:
                    raise e  # Raise the exception if it's not a conditional check failure

        print("CSV Content: ",csv_content)

        response_query = table.query(
            IndexName='ProcessingStatusIndex',
            KeyConditionExpression=boto3.dynamodb.conditions.Key('ProcessingStatus').eq('PENDING'),
        )

        QueuedEntries = response_query['Items']
        print("Items In Process: ", QueuedEntries)

        # Generate config.yaml and store in wgs-genomics-files/input
        config_data = generate_config(QueuedEntries)
        print("Config Data: ", config_data)

        s3_client.put_object(
            Bucket='wgs-genomics-yoyo458',
            Key='input/config.yaml',
            Body=yaml.dump(config_data).encode('utf-8')
        )

        response = initiate_snakemake_slurm()

        return {
            'statusCode': 200,
            'body': json.dumps('CSV processed successfully! HPC job submitted!')
        }
    
    except Exception as e:
        print(f"Error: {str(e)}")
        raise e
    
    

