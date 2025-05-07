module "api" {
    source = "aws-tf/parallelcluster/aws//modules/pcluster_api"
    api_version = var.api_version
    api_stack_name = var.api_stack_name
    parameters = {
        EnableIamAdminAccess = "true"
    }
}