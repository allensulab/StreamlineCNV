# StreamlineCNV

## About

All the pipeline tools are packaged into an OCI container image.
This image is available at [`ghcr.io/biomicrocenter/streamlinecnv:release`](https://github.com/orgs/BioMicroCenter/packages/container/package/streamlinecnv).
The `containers/` directory contains the Dockerfile used to create the image, alongside a folder containing the scripts which are included in the image.

To use this pipeline, you'll need to install Nextflow and either Singularity or Apptainer.

Install Singularity by following [the official installation instructions](https://docs.sylabs.io/guides/3.0/user-guide/installation.html) or install Apptainer by following [the official installation instructions](https://apptainer.org/docs/user/main/quick_start.html#installation).
Likewise, install Nextflow by following [the official installation instructions](https://www.nextflow.io/docs/latest/install.html#install-page).
If you are a user of a high-performance computing cluster (HPC), then you will need to check whether your cluster has Singularity/Apptainer and Nextflow installed already.
Please contact your system administrators for more information.

Once you've installed Singularity/Apptainer and Nextflow, clone this repository to your system: `git clone https://github.com/BioMicroCenter/StreamlineCNV && cd StreamlineCNV`.

This repository, including the OCI container image, is private. Thus, for either Singularity or Apptainer to pull the image you'll need to 1) have access to this repository, 2) create a personal access token for you account, and 3) configure Singularity or Apptainer to use your personal access token.

To create a personal access token, follow the [official GitHub account documentation](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#creating-a-personal-access-token-classic).

Then, add the following to whichever script you use (`run_local.sh`, `run_slurm.sh`, `run_container_no_nextflow.sh`):

```
# For Singularity:
export SINGULARITY_DOCKER_USERNAME="YOUR GITHUB USERNAME"
export SINGULARITY_DOCKER_PASSWORD="YOUR PERSONAL ACCESS TOKEN"  

# For Apptainer:
export APPTAINER_DOCKER_USERNAME="YOUR GITHUB USERNAME"
export APPTAINER_DOCKER_PASSWORD="YOUR PERSONAL ACCESS TOKEN"
```

Now, Singularity/Apptainer will use your GitHub credentials when attempting to pull the private StreamlineCNV OCI container image.

## Environment variables

Before you run the pipeline, you need to set the following environment variables to tell the pipeline where your data are located. 
You can do this in the provided example shell scripts.

| Variable        | Description                                                                                                           |
| :-------------- | :-------------------------------------------------------------------------------------------------------------------- |
| `DATA_DIR` | The folder where your reference files are located. This folder will be mounted at ‘/data’ inside the Singularity container. Defaults to ‘data/’. |
| `FASTQ_DATA` | The folder where your fastq files are located. Defaults to ‘data/fastq_files’. |
| `OUTDIR`        | The folder where you'd like results to be output. Defaults to `results`.                                                                     |
| `GENE_LIST`     | The gene list file. Defaults to `data/GeneListExample`.                                                                                  |
| `SAMPLE_INFO`   | The sample info file. Defaults to `data/SampleInfo`.                                                                                    |
| `CLUSTERING_LABEL`   | The file containing sample IDs to enter clustering analyses, as well as their labels such as cell types. Defaults to `data/ClusteringLabel`.                                                                                    |
| `QUEUE`         | (Optional) When using the `slurm` profile, the name of the Slurm queue you'll submit the job to. Defaults to `normal` |
| `CONTAINER_PROFILE`     | The container profile to use. Either `singularity` or `apptainer`. Defaults to `singularity` |
| `EXTRA_MOUNT`   | (Optional) Provide an extra mount point to the Singularity container. Defaults to `''`                                |

## Running the pipeline

### Run on a local workstation

Make sure you have Singularity/Apptainer and Nextflow installed.
Open the `run_local.sh` file and change the environment variables to reflect the location of the files on your computer.

Run `./run_local.sh`.
The results will be output to the directory you specified in `OUTDIR`. By default, this will be a folder named `results` in your current directory.

### Run on an HPC with Slurm

Make sure you have Singularity and Nextflow installed or loaded in to your environment.
Open the `run_slurm.sh` file and change the environment variables to reflect the location of the files on your server.

If the files are located on a storage array that your server mounts using NFS, set the `EXTRA_MOUNT` option to the path where the NFS mount starts.
For example, if your server mounts storage arrays at `/net/` and that's where your data are stored, set `EXTRA_MOUNT` to `/net:/net` to make those files available inside the container image at the same path that they would be on your server.

Submit the script to Slurm by running `sbatch run_slurm.sh`.

The results will be output to the directory you specified in `OUTDIR`. By default, this will be a folder named `results` in your current directory.

### Run without installing Nextflow, only Singularity/Apptainer

If for some reason you are able to install Singularity/Apptainer, but not Nextflow, you can still use this pipeline.
However, you will not be able to scatter your job across multiple compute nodes, limiting the functionality for those on high-performance computing environments.

Open the `run_container_no_nextflow.sh` file and change the environment variables to reflect the location of the files on your server.
Run `./run_container_no_nextflow.sh`.

Alternatively, you can submit the script to Slurm by running `sbatch run_container_no_nextflow.sh`.
Just make sure to edit the script to load Singularity/Apptainer into your HPC environment.
