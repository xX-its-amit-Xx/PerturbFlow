# AWS Batch deployment

PerturbFlow ships with an example AWS Batch job definition in
[`batch_job_definition.json`](batch_job_definition.json). It assumes:

- A Batch compute environment with at least one EC2 instance type providing
  16 vCPU / 64 GiB RAM (e.g. `r6i.4xlarge`).
- A scratch volume mounted at `/mnt/scratch` on the host (EBS or instance
  storage) — Perturb-seq matrices spill to disk for large screens.
- A container image published to your registry (the example points at GHCR).

## Register the job definition

```bash
# 1. Fill in your AWS account ID and image URI
sed -i "s/REPLACE_ACCOUNT_ID/$(aws sts get-caller-identity --query Account --output text)/g" \
  aws/batch_job_definition.json

# 2. Register
aws batch register-job-definition \
  --cli-input-json file://aws/batch_job_definition.json
```

## Submit a run

```bash
aws batch submit-job \
  --job-name perturbflow-replogle-2022 \
  --job-queue perturbflow-queue \
  --job-definition perturbflow \
  --parameters config_path=/scratch/run/config.yaml,output_path=/scratch/run/perturbflow
```

## IAM

Two roles referenced in the job definition:

- `PerturbFlowExecutionRole` — pulls the container image, writes logs to
  CloudWatch. Standard `AmazonECSTaskExecutionRolePolicy`.
- `PerturbFlowJobRole` — what the running pipeline assumes. Grant
  `s3:GetObject`/`s3:PutObject` on the bucket holding your input matrices
  and output artifacts.

Do **not** reuse the execution role as the job role — the pipeline should
only see the data buckets it needs, not the registry pull credentials.
