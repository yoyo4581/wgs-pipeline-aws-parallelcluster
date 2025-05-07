import json
import urllib.parse
import boto3
import os
import time
import csv
import io
import yaml


def lambda_handler(event, context):
    print("Received event: " + json.dumps(event, indent=2))

    user = event['user']
    account = event['slurm_account']
    partition = 'demo'

    

