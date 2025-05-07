module "required_infra" {
  source = "aws-tf/parallelcluster/aws//modules/required_infra"

  prefix               = var.prefix
  vpc_cidr_block       = var.vpc_cidr_block
  public_subnet_cidrs  = var.public_subnet_cidrs
  private_subnet_cidrs = var.private_subnet_cidrs

  public_subnet_az  = var.public_subnet_az
  private_subnet_az = var.private_subnet_az
}