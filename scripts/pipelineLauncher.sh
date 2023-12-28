#!/bin/bash
eval "$(conda shell.bash hook)"
conda activate ngs-pipeline-launcher
RUN=$1
EMAIL=${2-None}
sbatch <<EOT
#!/bin/bash

#SBATCH --job-name=PipelineLauncher
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00:00
#SBATCH --mem=1G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$EMAIL
#SBATCH --output=$RUN"%x_out.txt"
#SBATCH --error=$RUN"%x_error.txt"
#SBATCH --partition=vm-cpu

python /nfs/APL_Genomics/apps/production/ngs-pipeline-launcher/scripts/pipelineLauncher.py -r $RUN -e $EMAIL

exit 0
EOT 