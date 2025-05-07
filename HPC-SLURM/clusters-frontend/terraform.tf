terraform {
  required_providers {
    aws = {
      source  = "hashicorp/aws"
      version = "~> 5.82.0"
    }
    awscc = {
      source  = "hashicorp/awscc"
      version = ">= 1.34.0"
    }
    aws-parallelcluster = {
      source  = "aws-tf/aws-parallelcluster"
      version = "1.1.0" 
    }
  }
}