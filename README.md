# seq-sample-split

Overall workflow:

1. Start the sequencing run
2. Create the sample sheet using the template
3. Start the script and it'll check every 15 minutes to see if sequencing is done
4. When done, the sequencing output files are moved to specific directories based on the sample sheet
5. A SLURM file is created for each directory/pipeline and sent to the scheduler