
variable "region" {
  description = "The region used to deploy the infrastructure."
  type        = string
  default     = "us-east-1"
}

variable "profile" {
  type        = string
  description = "The AWS profile used to deploy the infrastructure."
  default     = null
}

## GENERAL ###############

variable "prefix" {
  type        = string
  description = "The prefix used for resource names."
  default     = "pcluster"
}

## VPCs ##################

variable "vpc_cidr_block" {
  type        = string
  description = "The CIDR block of the VPC the cluster nodes will be created in. The public and private subnet CIDR blocks should fall within this block."
  default     = "10.0.0.0/16"
}

## SUBNETs ###############

variable "public_subnet_cidrs" {
  type        = list(any)
  description = "The list of CIDR blocks to be used for public subnets. Has to be in the VPC CIDR block. Cannot conflict with private subnets."
  default     = ["10.0.1.0/24"]
}

variable "private_subnet_cidrs" {
  type        = list(any)
  description = "The list of CIDR blocks to be used for private subnets. Has to be in the VPC CIDR block. Cannot conflict with public subnets."
  default     = ["10.0.2.0/24"]
}

variable "public_subnet_az" {
  type        = string
  description = "The availability zone to create the public subnets in."
  default     = "us-east-1a"
}

variable "private_subnet_az" {
  type        = string
  description = "The availability zone to create the private subnets in."
  default     = "us-east-1b"
}