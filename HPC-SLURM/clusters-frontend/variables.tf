/**
 *  Copyright 2024 Amazon.com, Inc. or its affiliates. All Rights Reserved.
 *
 *  Licensed under the Apache License, Version 2.0 (the "License"). You may not
 *  use this file except in compliance with the License. A copy of the License is
 *  located at
 *
 *  http://aws.amazon.com/apache2.0/
 *
 *  or in the "LICENSE.txt" file accompanying this file. This file is distributed
 *  on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, express or
 *  implied. See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 */

variable "region" {
  description = "The region to deploy the ParallelCluster API."
  type        = string
  default     = "us-east-1"
}

variable "profile" {
  type        = string
  description = "The AWS profile used to deploy resources."
  default     = null
}

variable "public_subnet_az" {
  type        = string
  description = "The availability zone to create the public subnets in."
  default     = "us-east-1a"
}

variable "private_subnet_az" {
  type        = string
  description = "The availability zone to create the private subnets in."
  default     = "us-east-1a"
}

variable "api_version" {
  type        = string
  description = "The version of the ParallelCluster API to deploy."
  default = "3.13.0"
}

variable "api_stack_name" {
  type        = string
  description = "The name of the CloudFormation stack use to deploy the ParallelCluster API."
  default     = "ParallelCluster"
}

variable "project_name" {
  description = "Name of the project to be used in resource naming"
  type        = string
  default = "wgs-hpc"

  validation {
    condition     = can(regex("^[a-zA-Z][a-zA-Z0-9-]*$", var.project_name))
    error_message = "Project name must start with a letter and can only contain letters, numbers, and hyphens."
  }
}


variable "db_name" {
  type        = string
  description = "The accounting database name."
  default     = "accountingdb"
}

variable "db_username" {
  description = "Slurm database username"
  type        = string
  default    = "yoyo458"
}