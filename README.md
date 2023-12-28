# ngs-pipeline-launcher
 [![Lifecycle: WIP](https://img.shields.io/badge/lifecycle-WIP-yellow.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) [![Contributions welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg?style=flat)](https://github.com/provlab-bioinfo/ngs-pipeline-launcher/issues) [![License: GPL3](https://img.shields.io/badge/license-GPL3-lightgrey.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html) [![minimal Python version: 3.0](https://img.shields.io/badge/Python-3.10-6666ff.svg)](https://www.python.org/) [![Package Version = 0.0.1](https://img.shields.io/badge/Package%20version-0.0.1-orange.svg?style=flat-square)](https://github.com/provlab-bioinfo/ngs-pipeline-launcher/blob/main/NEWS) [![Last-changedate](https://img.shields.io/badge/last%20change-2023--10--31-yellowgreen.svg)](https://github.com/provlab-bioinfo/ngs-pipeline-launcher/blob/main/NEWS)

A launcher tool for sorting a sequencing run by pathogen and automatic initialization of their downstream pipelines. Build for routine NGS at [ProvLab](https://github.com/provlab-bioinfo/pathogenseq).

## Table of Contents

- [Quick-Start Guide](#quick-start%guide)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Input](#input)
- [Output](#output)
- [Pipeline Usage](#pipeline-usage)
- [References](#references)

## Quick-Start Guide

When ```PipelineWorksheet.xlsx``` is filled in, launch with:

```bash
conda activate ngs-pipeline-launcher
run_pipeline_launcher [/path/to/run/] [E-mail]
```

## Dependencies

[Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) is required to build the [environment](/environments/environment.yml) with the necessary workflow dependencies. To create the environment and setup the run command:
```bash
conda env create -f ./environments/environment.yml
conda activate ngs-pipeline-launcher
cd $CONDA_PREFIX
mkdir -p ./etc/conda/activate.d
mkdir -p ./etc/conda/deactivate.d
touch ./etc/conda/activate.d/env_vars.sh
touch ./etc/conda/deactivate.d/env_vars.sh
printf '#!/bin/sh\nrun_pipeline_launcher() { bash /path/to/pipelineLauncher.sh $1 $2; }\nexport -f run_pipeline_launcher' > ./etc/conda/activate.d/env_vars.sh
printf '#!/bin/sh\nunset run_pipeline_launcher' > ./etc/conda/deactivate.d/env_vars.sh
```

## Input

#### PipelineWorksheet.xlsx:

This file species the run information in ```[Header]```, the output directories in ```[Directories]```, the pipelines to use in ```[Pipelines]```, and the sample information in ```[Samples]```.

## Output

The output files are determined by the particular pipelines being run, and will be located in the folders specified by ```[Directories]``` in the ```PipelineWorksheet.xlsx```

## More information

Please see ```\Laboratory Documents\Cross-Site\Research\Starting Pathogen Analysis via the Pipeline Launcher``` on Paradigm