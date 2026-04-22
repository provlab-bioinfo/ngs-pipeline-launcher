# This will create the environment and the system command to run the launcher
conda env create -f ./environments/environment.yml -y
conda activate ngs-pipeline-launcher
cd $CONDA_PREFIX
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
touch ./etc/conda/deactivate.d/env_vars.sh
printf '#!/bin/sh\nrun_pipeline_launcher() { bash /path/to/pipelineLauncher.sh $1 $2; }\nexport -f run_pipeline_launcher\nmove_sequencing_files() { python /path/to/moveToAPLGenomics.py --p $1 --d $2; }\nexport -f move_sequencing_files' >> ./etc/conda/activate.d/env_vars.sh
printf '#!/bin/sh\nunset run_pipeline_launcher\unset move_sequencing_files' >> ./etc/conda/deactivate.d/env_vars.sh