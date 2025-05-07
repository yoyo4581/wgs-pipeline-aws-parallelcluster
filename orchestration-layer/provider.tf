provider "aws-parallelcluster" {
  region = "us-east-1"
}

provider "aws" {
    region = "us-east-1"
    shared_credentials_files = ["/Users/Yahya/.aws/credentials"]
}
