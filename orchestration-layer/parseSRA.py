import boto3
import json

def get_pending_samples(table_name):
    dynamodb = boto3.resource("dynamodb")
    table = dynamodb.Table(table_name)
    response = table.scan(
        FilterExpression="attribute_not_exists(processed) OR processed = :val",
        ExpressionAttributeValues={":val": False}
    )
    return response['Items']

if __name__ == "__main__":
    samples = get_pending_samples("SRATable")
    with open("pending_samples.json", "w") as f:
        json.dump(samples, f)
