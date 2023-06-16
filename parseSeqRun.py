import pandas as pd, os, searchTools as st, shutil, io, time, fileinput, subprocess
from configparser import ConfigParser
pd.options.mode.chained_assignment = None  # default='warn'
os.chdir(os.path.dirname(__file__))

def getSampleSheetDataVars(path:str, section:str):
    """Generates a dictionary from the first two columns of a [HEADER] section
    :param path: The path to the reference file
    :param section: The name of the section
    :return: A dictionary containing variable keys
    """    
    cfg = ConfigParser(allow_no_value=True)
    cfg.optionxform = str
    cfg.read(path)
    dict = {k:v for k, *v in map(lambda x: str.split(x,sep=","), cfg[section])}
    dict = {k:v[0] for k, v in dict.items() if v} # Removes blank keys
    return (dict)

def getSampleSheetDataFrame(path:str, section:str):
    """Generates a DataFrame from a [HEADER] section
    :param path: The path to the reference file
    :param section: The name of the section
    :return: A DataFrame representing the sections
    """ 
    cfg = ConfigParser(allow_no_value=True)
    cfg.optionxform = str
    cfg.read(path)
    buf = io.StringIO()
    buf.writelines('\n'.join(row.rstrip(',') for row in cfg[section]))
    buf.seek(0)
    return (pd.read_csv(buf))

def isRunCompleted(path:str, seqType: str):
    """Checks whether a sequencing run is completed. 
    For Illumina, checks for the 'CompletedJobInfo.xml' file. 
    For Nanopore, checks for the 'final_summary_*.txt' file.
    :param path: The path to the experiment directory
    :param seqType: The type of sequencing. Either 'Illumina' or 'Nanopore'
    :return: True if complete, False if not
    """    
    file = ""
    if seqType == "Nanopore":
        file = "final_summary_*.txt"
    elif seqType == "Illumina":
        file = "CompletedJobInfo.xml"
    return False if not len(st.findFile(os.path.join(path,"**",file))) else True

def generateSLURM(SLURM:str, jobName: str, outputDir: str, command: str):
    """Generates a SLURM command file based on a template
    :param SLURM: Path to the template SLURM file
    :param jobName: The name of the job
    :param outputDir: The directory for output/error files
    :param command: The path to the new SLURM file
    """
    file = open(SLURM, "rt")
    data = file.read()
    file.close()
    data = data.replace("[JOB_NAME]", jobName)
    data = data.replace("[OUTPUT_DIR]", os.path.join(outputDir,"SLURM_out.txt"))
    data = data.replace("[ERROR_DIR]", os.path.join(outputDir,"SLURM_error.txt"))
    data = data.replace("[RUN_DIR]", outputDir)
    outFile = os.path.join(outputDir,os.path.basename(SLURM))
    file = open(outFile, "wt+")
    file.write(data)
    file.write("\n\n"+command)
    file.close()
    return(outFile)

sampleSheetPath = "SampleSheet.csv"
allSamples = getSampleSheetDataFrame(sampleSheetPath, "Samples")
directories = getSampleSheetDataVars(sampleSheetPath, "Directories")
pipelines = getSampleSheetDataVars(sampleSheetPath, "Pipelines")
header = getSampleSheetDataVars(sampleSheetPath, "Header")
basePath = header["Run_Dir"]
SLURM = "SLURM.batch"

# Check if run is finished sequencing
while not isRunCompleted(basePath, header["Seq_Type"]):
    print("waiting...")
    time.sleep(15*60)

# # Split sequencing run into respective folders
# for group in set(allSamples["Sample_Group"].values):
#     if group == "Ignore": continue
#     outDir = directories[group]
#     print("Moving " + group + " to " + outDir)
#     excludeSamples = allSamples.loc[allSamples['Sample_Group'] != group]["Sample_ID"].values.tolist()
#     excludeSamples = ",".join(excludeSamples).split(",")
#     excludeSamples = ["*"+sample+"*" for sample in excludeSamples]
#     excludeSamples = excludeSamples + ["*fail*","*skip*","*unclassified*"]
#     shutil.copytree(basePath, outDir, ignore=shutil.ignore_patterns(*excludeSamples))

# Setup pipeline
for group in set(allSamples["Sample_Group"].values):
    if group == "Ignore": continue
    # Parse controls
    samples = allSamples.loc[allSamples['Sample_Group'] == group]
    ctrls = samples.dropna()
    samples = samples[samples['Control'].isna()]

    negCtrls = ctrls.loc[ctrls['Control'] == "negative"]
    if (len(negCtrls)):
        negCtrls = ",".join(map(str,negCtrls["Sample_Pos"].values.tolist()))

    posCtrls = ctrls.loc[ctrls['Control'] != "negative"]
    if (len(posCtrls)):
        posCtrls["Control"] = posCtrls["Sample_Pos"].astype(str) +","+ posCtrls["Control"].astype(str)
        posCtrls = " ".join(posCtrls["Control"].values.tolist())

    # Create the SLURM command

    symlink = lambda dir,link: "ln -s {} {}".format(st.findFile(os.path.join(directories[group],"**",dir))[0],)
    symlink = lambda dir,link: "ln -s {} {}".format(st.findFile(os.path.join("./**",dir))[0],os.path.join(directories[group],link))

    def symlink(dir,link): # Creates symlink command
        file = st.findFile(os.path.join(directories[group],"**",dir))[0]
        file = os.path.relpath(file,directories[group])
        link = os.path.join(directories[group],link)
        link = os.path.relpath(link,directories[group])
        return "ln -s {} {}".format(file,link)

    commands = []
    if (group == "ncov"):
        # Do symlinks
        if header["Seq_Type"] == "Nanopore":
            commands.append(symlink("fast5_pass","fast5"))   
            commands.append(symlink("fastq_pass","gup_out"))   
        elif header["Seq_Type"] == "Illumina":
            commands.append(symlink("fastq","Fastq")) 

        # Go to parent dir
        parentDir = os.path.dirname(directories[group].rstrip("/")) + "/"
        commands.append("\ncd {}\n".format(parentDir))

        # Run pipeline
        command = "python {} -d {} -r {} -b 2".format(pipelines[group],parentDir,header["Run_Name"])
        if len(posCtrls): command = command + " -p {}".format(posCtrls)
        if len(negCtrls): command = command + " -c {}".format(negCtrls)
        commands.append(command)

    elif (group == "ncov-ww"):
        # Do symlinks
        if header["Seq_Type"] == "Nanopore":
            commands.append(symlink("fast5_pass","fast5"))   
            commands.append(symlink("fastq_pass","gup_out"))   
        elif header["Seq_Type"] == "Illumina":
            commands.append(symlink("fastq","Fastq")) 

        # Go to parent dir
        parentDir = os.path.dirname(directories[group].rstrip("/")) + "/"
        commands.append("\ncd {}\n".format(parentDir))

        # Run pipeline
        command = "python {} -d {} -r {} -b 2 -f".format(pipelines[group],parentDir,header["Run_Name"])
        if len(posCtrls): command = command + " -p {}".format(posCtrls)
        if len(negCtrls): command = command + " -c {}".format(negCtrls)
        commands.append(command)

    # Generate the SLURM file
    SLURMfile = generateSLURM(SLURM, header["Run_Name"]+"_"+group, directories[group], "\n".join(commands))
    subprocess.run(["sbatch",SLURMfile,"-v"])
