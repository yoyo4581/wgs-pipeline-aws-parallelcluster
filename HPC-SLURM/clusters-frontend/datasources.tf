data "terraform_remote_state" "vpc" {
    backend = "local"
    config = {
        path = "../clusters-backend/infrastructure/terraform.tfstate"
    }
}

data "terraform_remote_state" "api" {
    backend = "local"
    config = {
        path = "../clusters-backend/api/terraform.tfstate"
    }
}

data "aws_lambda_function" "pcluster_function" {
  function_name = local.lambda_function_name
}

data "aws_iam_role" "api_user_role" {
  name = split("/", data.terraform_remote_state.api.outputs.pcluster_api_stack_outputs.ParallelClusterApiUserRole)[1]
}


data "aws_secretsmanager_secret" "slurm_db_secret" {
  name = "slurmdb-credentials"
}

data "aws_secretsmanager_secret_version" "db_password" {
  secret_id = data.aws_secretsmanager_secret.slurm_db_secret.id
}


data "aws_security_group" "default_sg" {
  vpc_id = data.terraform_remote_state.vpc.outputs.vpc_id
  name = "default"
}
