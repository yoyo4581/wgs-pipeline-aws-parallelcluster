provider "aws" {
  region = var.region
  shared_credentials_files = ["/Users/Yahya/.aws/credentials"]
}

provider "awscc" {
  region = var.region
}

provider "aws-parallelcluster" {
  region   = var.region
  endpoint = data.terraform_remote_state.api.outputs.pcluster_api_stack_outputs.ParallelClusterApiInvokeUrl
  role_arn = data.terraform_remote_state.api.outputs.pcluster_api_stack_outputs.ParallelClusterApiUserRole
}


