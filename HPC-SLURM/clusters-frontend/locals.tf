locals {
  config_vars = {
    subnet   = data.terraform_remote_state.vpc.outputs.public_subnets_ids[0],
    key_pair = data.terraform_remote_state.vpc.outputs.key_pair_id,
    region   = var.region
    db_uri = aws_db_instance.accountDB.endpoint
    db_username = aws_db_instance.accountDB.username
    secret_arn = data.aws_secretsmanager_secret.slurm_db_secret.arn
    db_name = aws_db_instance.accountDB.db_name
    profile = var.profile
    project_name = var.project_name
    default_sg = data.aws_security_group.default_sg.id
    ami_id = "ami-09dd073269896fe81"
  }
}

locals {
  lambda_function_name = split(":", data.terraform_remote_state.api.outputs.pcluster_api_stack_outputs.ParallelClusterLambdaArn)[6]
}


locals {
  lambda_role_arn = data.aws_lambda_function.pcluster_function.role
  lambda_role_name = split("/", local.lambda_role_arn)[1]
}