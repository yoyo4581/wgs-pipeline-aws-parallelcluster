module "pcluster" {
  source                = "aws-tf/parallelcluster/aws"
  deploy_pcluster_api   = false
  deploy_required_infra = false
  region                = var.region
  api_stack_name        = var.api_stack_name
  api_version           = var.api_version
  template_vars         = local.config_vars
  config_path           = "files/pcluster-example-config.yaml"
}




resource "aws_db_instance" "accountDB" {
  allocated_storage    = 5
  db_name              = "accountdb"
  engine               = "mysql"
  engine_version       = "8.0"
  instance_class       = "db.t3.micro"
  username             = var.db_username
  password             = data.aws_secretsmanager_secret_version.db_password.secret_string
  parameter_group_name = "default.mysql8.0"
  skip_final_snapshot  = true
  vpc_security_group_ids = [aws_security_group.rds_sg.id]
  db_subnet_group_name = aws_db_subnet_group.rds_subnet_group.name
}

resource "aws_security_group" "rds_sg" {
  name        = "rds-security-group"
  description = "Security group for RDS instance"
  vpc_id      = data.terraform_remote_state.vpc.outputs.vpc_id

  ingress {
    from_port       = 3306
    to_port         = 3306
    protocol        = "tcp"
    security_groups = [data.aws_security_group.default_sg.id]
  }
}

resource "aws_db_subnet_group" "rds_subnet_group" {
  name       = "accountdb-subnet-group"
  subnet_ids = flatten([
    data.terraform_remote_state.vpc.outputs.public_subnets_ids,
    data.terraform_remote_state.vpc.outputs.private_subnets_ids
  ])
}


resource "aws_iam_role_policy_attachment" "s3_fullaccess" {
  role       = data.aws_iam_role.api_user_role.name
  policy_arn = "arn:aws:iam::aws:policy/AmazonS3FullAccess"
}
resource "aws_iam_role_policy_attachment" "secrets_access" {
  role       = data.aws_iam_role.api_user_role.name
  policy_arn = "arn:aws:iam::aws:policy/SecretsManagerReadWrite"
}


resource "aws_iam_role_policy_attachment" "s3_fullaccess_lambda" {
  role       = local.lambda_role_name
  policy_arn = "arn:aws:iam::aws:policy/AmazonS3FullAccess"
}
resource "aws_iam_role_policy_attachment" "secrets_access_lambda" {
  role       = local.lambda_role_name
  policy_arn = "arn:aws:iam::aws:policy/SecretsManagerReadWrite"
}

resource "aws_iam_role_policy" "allow_attach_role_policy" {
  name = "AllowAttachRolePolicy"
  role = local.lambda_role_name

  policy = jsonencode({
    Version = "2012-10-17",
    Statement = [
      {
        Effect = "Allow",
        Action = [
          "iam:AttachRolePolicy",
          "iam:DetachRolePolicy",
        ],
        Resource = "*"
      }
    ]
  })
}
