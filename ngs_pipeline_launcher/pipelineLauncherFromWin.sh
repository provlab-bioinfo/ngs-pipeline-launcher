#!/bin/bash

file=$1
source activate ngs-pipeline-launcher
run_pipeline_launcher $file

sq() {
	squeue -o "%.20i %.2t | %.40u %.40j | %.10M %.5D %.4C %.7m | %.14R" "$@" \
	| sed \
		-e "s/ MIN_MEM / MEM /g" \
		-e "s/ JOBID / JOB_ID /g" \
		-e "s/ NAME / JOB_NAME /g" \
	| column -t \
	| sed \
	-e "s/ R / $(echo -e '\033[92mR\033[0m') /g" \
	-e "s/ PD / $(echo -e '\033[93mPD\033[0m') /g" \
	-e "s/|/$(tput bold; tput setaf 233)|$(tput sgr0)/g" \
	-e "1s/\x1B\[[0-9;]*[a-zA-Z]//g" \
	-e "1s/^/$(tput bold; tput setab 243;)/;1s/$/\x1B[0m/" \
	-e "s/^/        /" \
	-e "1s/^/        🖥️  SLURM Status as of $(date +"%H:%M:%S %p") on $(date +"%b %d")\n/" \
	-e "1s/^/\n/" \
	-e '$s/$/\n/' \
	-e "s/\b$USER\b/$(tput bold; tput setaf 7)$USER$(tput sgr0)/g" \
	-e "s/ NODELIST(REASON)/ NODE          /g" 
}

cwatch() {     
	local interval=1   
	if [[ "$1" == "-n" ]]
	then interval=$2       
	shift 2  
	fi
	while true
	do refreshed=$(eval "$@")
	clear
	#echo "Every ${interval}s: $@"       
	#echo "----------------------"      
	echo "$refreshed"
	sleep "$interval"  
	done
}

cwatch sq