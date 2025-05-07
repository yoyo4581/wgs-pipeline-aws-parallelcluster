#GEOAccession, SampleName, SRAStudy
resource "aws_dynamodb_table" "samples-table" {
  name           = "SRATable"
  billing_mode   = "PAY_PER_REQUEST"
  hash_key       = "ExperimentId"
  range_key      = "SRRId"

  attribute {
    name = "SRRId"
    type = "S"
  }

  attribute {
    name = "ExperimentId"
    type = "S"
  }

  attribute {
    name = "CellLine"
    type = "S"
  }
  attribute {
    name = "ProcessingStatus"
    type = "S"
  }
  # attribute {
  #   name = "MetadataBlob"
  #   type = "S"
  # }

  # attribute {
  #   name = "AvgSpotLength"
  #   type = "N"
  # }

  # attribute {
  #   name = "Platform"
  #   type = "S"
  # }
  # attribute {
  #   name = "LibraryLayout"
  #   type = "S"
  # }

  ttl {
    attribute_name = "TimeToExist"
    enabled        = true
  }

  global_secondary_index {
    name               = "CellLineIndex"
    hash_key           = "CellLine"
    range_key          = "ExperimentId"
    projection_type    = "ALL"
  }

  global_secondary_index {
    name               = "ProcessingStatusIndex"
    hash_key           = "ProcessingStatus"
    range_key          = "ExperimentId"
    projection_type    = "ALL"
  }

  tags = {
    Name        = "SRATable"
    Environment = "production-wgs-hpc"
  }
}


#--------------------------------------------------------------------

resource "aws_s3_bucket" "snakemake-folder" {
    bucket = "wgs-snakemake-files-yoyo458"
}

# Upload the Snakemake file
resource "aws_s3_object" "snakemake_file" {
  bucket = aws_s3_bucket.snakemake-folder.bucket
  key    = "Snakefile"
  source = "./buckets/snakemake-files/Snakefile"
  acl    = "private"
  etag   = filemd5("./buckets/snakemake-files/Snakefile")
}

resource "aws_s3_object" "snakemake-additional" {
  for_each = fileset("./buckets/snakemake-files/snakemake_additional_files", "*.smk")

  bucket = aws_s3_bucket.snakemake-folder.bucket
  key    = "snakemake_additional_files/${each.value}"
  source = "./buckets/snakemake-files/snakemake_additional_files/${each.value}"
  # etag makes the file update when it changes; see https://stackoverflow.com/questions/56107258/terraform-upload-file-to-s3-on-every-apply
  etag   = filemd5("./buckets/snakemake-files/snakemake_additional_files/${each.value}")
}

resource "aws_s3_object" "snakemake-setup" {
  bucket = aws_s3_bucket.snakemake-folder.bucket
  key = "ami-build/setup_pipeline.sh"
  source = "./buckets/setup-install/setup_pipeline.sh"
  acl = "private"
  etag = filemd5("./buckets/setup-install/setup_pipeline.sh")
}

#----------------------------------------------------------------------
data "aws_s3_bucket" "wgs-genomics-bucket" {
    bucket = "wgs-genomics-yoyo458"
}



resource "aws_s3_bucket" "input-SRA" {
    bucket = "input-sra-yoyo458"
}



resource "aws_s3_object" "input-SRA-files" {
  for_each = fileset("./SRAInput/input", "*.csv")

  bucket = aws_s3_bucket.input-SRA.bucket
  key    = "input/${each.value}"
  source = "./SRAInput/input/${each.value}"
  # etag makes the file update when it changes; see https://stackoverflow.com/questions/56107258/terraform-upload-file-to-s3-on-every-apply
  etag   = filemd5("./SRAInput/input/${each.value}")
}


resource "aws_s3_object" "input-files" {
  bucket = data.aws_s3_bucket.wgs-genomics-bucket.bucket
  key    = "reference/.keep"  # Creating a dummy file to represent the "temp" prefix
  acl    = "private"
  content = ""  # Empty content
}
resource "aws_s3_object" "reference-files" {
  bucket = data.aws_s3_bucket.wgs-genomics-bucket.bucket
  key    = "reference/.keep"  # Creating a dummy file to represent the "temp" prefix
  acl    = "private"
  content = ""  # Empty content
}
resource "aws_s3_object" "temp-files" {
  bucket = data.aws_s3_bucket.wgs-genomics-bucket.bucket
  key    = "temp/.keep"  # Creating a dummy file to represent the "temp" prefix
  acl    = "private"
  content = ""  # Empty content
}
resource "aws_s3_object" "output-files" {
  bucket = data.aws_s3_bucket.wgs-genomics-bucket.bucket
  key    = "output/.keep"  # Creating a dummy file to represent the "temp" prefix
  acl    = "private"
  content = ""  # Empty content
}

data "aws_caller_identity" "current" {}

output "account_id" {
  value = data.aws_caller_identity.current.account_id
  sensitive = true
}

data "archive_file" "process-s3-sra" {
    type = "zip"
    source_dir = "./lambda_functions/ProcessSRA"
    output_path = "./lambda_functions/ProcessSRA/process_s3_sra.zip"
    excludes = ["*.zip"]
}

data "archive_file" "process-s3-sra_layer" {
    type = "zip"
    source_dir = "./layers/lambda_layers/"
    output_path = "./layers/lambda_layer_yaml.zip"
}

resource "aws_iam_policy" "LambdaDynamoSRAPolicy" {
  name = "LambdaDynamoSRAPolicy"
  path = "/"
  description = "Policy for Lambda function to access DynamoDB and S3"

  policy = jsonencode({
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "dynamodb:PutItem",
                "dynamodb:UpdateItem",
                "dynamodb:GetItem",
                "dynamodb:Query",
                "dynamodb:Scan"
            ],
            "Resource": "arn:aws:dynamodb:us-east-1:302263065877:table/SRATable"
        },
        {
          "Effect": "Allow",
          "Action": [
            "dynamodb:Query",
            "dynamodb:Scan"
          ],
          "Resource": "arn:aws:dynamodb:us-east-1:302263065877:table/SRATable/*" 
        },
        {
            "Effect": "Allow",
            "Action": [
                "logs:CreateLogGroup",
                "logs:CreateLogStream",
                "logs:PutLogEvents"
            ],
            "Resource": "*"
        }
    ]
  })
}

resource "aws_iam_role_policy_attachment" "lambda_dynamodb_policy_attachment" {
  role = "WGS-Trigger-LambdaRole"
  policy_arn = aws_iam_policy.LambdaDynamoSRAPolicy.arn
}

resource "aws_lambda_layer_version" "pyyaml_layer" {
  filename = data.archive_file.process-s3-sra_layer.output_path
  layer_name = "yaml_layer"

  compatible_runtimes = ["python3.12"]
}

resource "aws_lambda_function" "process-s3-sra" {
    description = "Lambda function to process S3 SRA files into DynamoDB and trigger Snakemake via Fargate task"

    function_name = "process-s3-sra"
    role = "arn:aws:iam::${data.aws_caller_identity.current.account_id}:role/WGS-Trigger-LambdaRole"
    handler = "process_s3_sra.lambda_handler"
    runtime = "python3.12"
    layers= [aws_lambda_layer_version.pyyaml_layer.arn]
    filename = data.archive_file.process-s3-sra.output_path
    source_code_hash = data.archive_file.process-s3-sra.output_base64sha256
    tags = {
        "lambda-console:blueprint" = "s3-get-object-python"
    }
}

resource "aws_lambda_permission" "allow_bucket" {
    action = "lambda:InvokeFunction"
    function_name = aws_lambda_function.process-s3-sra.id
    source_account = "${data.aws_caller_identity.current.account_id}"
    principal = "s3.amazonaws.com"
    source_arn = aws_s3_bucket.input-SRA.arn
}

resource "aws_s3_bucket_notification" "upload_sra_notification" {
  bucket = aws_s3_bucket.input-SRA.id
  eventbridge = true

  lambda_function {
    lambda_function_arn = aws_lambda_function.process-s3-sra.arn
    events              = ["s3:ObjectCreated:Put"]

    # Optional filters:
    filter_prefix = "input/"
    filter_suffix = ".csv"
  }

  depends_on = [aws_lambda_permission.allow_bucket]
}



