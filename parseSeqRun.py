import pandas as pd, os, searchTools as st, shutil, io, time, fileinput, subprocess, argparse, tempfile, pathlib
from configparser import ConfigParser
from datetime import datetime
pd.options.mode.chained_assignment = None  # default='warn'
os.chdir(os.path.dirname(__file__))

defaultSampleSheet = "/nfs/Genomics_DEV/projects/alindsay/Projects/seq-sample-split/PipelineWorksheet.xlsx"
validPipelines = ["ncov","ncov-ww"]

barcodeCol = "Barcode"
samplePosCol = "Sample_Pos"

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
    dict = {k:v[0] for k, v in dict.items() if v} # Removes blank keys and keep only first column after var
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
    df = pd.read_csv(buf)
    return (df)

def isRunCompleted(path:str, seqType: str):
    """Checks whether a sequencing run is completed. 
    For Illumina, checks for the 'CompletedJobInfo.xml' file. 
    For Nanopore, checks for the 'final_summary_*.txt' file.
    :param path: The path to the experiment directory
    :param seqType: The type of sequencing. Either 'Illumina' or 'Nanopore'
    :return: True if complete, False if not
    """    
    if not os.path.exists(path):
        raise Exception("Run directory does not exist")

    file = ""
    if seqType.lower() == "nanopore":
        file = "final_summary_*.txt"
    elif seqType.lower() == "illumina":
        file = "CompletedJobInfo.xml"

    found = st.findFile(os.path.join(path,"**",file))

    if (len(found)):
        print(f"   Found {file} at {found}")
        return True
    else:
        return False

    # return False if not len(st.findFile(os.path.join(path,"**",file))) else True

def generateSLURM(SLURM:str, jobName: str, runName: str, outputDir: str, command: str, email: str = None):
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
    data = data.replace("[OUTPUT_DIR]", os.path.join(outputDir,runName+"_out.txt"))
    data = data.replace("[ERROR_DIR]", os.path.join(outputDir,runName+"_error.txt"))
    data = data.replace("[EMAIL]", email)
    data = data.replace("[MAIL_TYPE]", "NONE" if email is None else "ALL")      
    data = data.replace("[RUN_DIR]", outputDir)
    outFile = os.path.join(outputDir,runName+"_SLURM.batch")
    file = open(outFile, "wt+")
    file.write(data)
    file.write("\n\n"+command)
    file.close()
    return(outFile)

print(f"Pipeline launcher started at: {datetime.now().strftime('%H:%M:%S')}")

# Import the sample sheet
parser = argparse.ArgumentParser(description='APL NGS Pipeline Launcher')
parser.add_argument("-w", "--worksheet", help="Path to the pipeline worksheet.", default = defaultSampleSheet)
parser.add_argument("-e", "--email", help="Notify status alerts by e-mail.", default = None)
args = parser.parse_args()
sampleSheetPath = args.worksheet
email = args.email

with tempfile.NamedTemporaryFile() as sampleSheet:
    if pathlib.Path(args.worksheet).suffix == ".xlsx":
        df = pd.read_excel(args.worksheet)
        sampleSheetPath = sampleSheet.name
        df.to_csv(sampleSheetPath, index=False)
    else: sampleSheetPath = args.worksheet

    allSamples = getSampleSheetDataFrame(sampleSheetPath, "Samples")
    directories = getSampleSheetDataVars(sampleSheetPath, "Directories")
    pipelines = getSampleSheetDataVars(sampleSheetPath, "Pipelines")
    header = getSampleSheetDataVars(sampleSheetPath, "Header")

basePath = header["Run_Dir"]
SLURM = "SLURM_template.batch"

# Check for appropriate inputs
if header["Seq_Type"].lower() not in ["nanopore","illumina"]:
    raise Exception(f"Seq_Type must be either 'Nanopore' or 'Illumina' in the pipeline worksheet. Found: {header['Seq_Type']}")

if header["Base_Call"].lower() not in ["yes","no"]:
    raise Exception(f"Base_Call must be either 'Yes' or 'No' in the pipeline worksheet. Found: {header['Base_Call']}")

groups = sorted(set(allSamples["Sample_Group"].dropna().values))
for group in groups:

    if group == "ignore":
        continue

    if group not in validPipelines:
        raise Exception(f"Pipeline '{group}' is not a valid option")
    
    if not os.path.exists(pipelines[group]):
        if (pipelines[group] != "ignore"):
            raise Exception(f"Pipeline for '{group}' does not exist at '{pipelines[group]}'")

    if os.path.exists(directories[group]):
        if (directories[group] != "ignore"):
            if (len(directories[group]) != 0):
                raise Exception(f"Directory for '{group}' at '{directories[group]}' already exists and is not empty. Please choose empty or non-existing directory.")

print("Checking for sequencing completion file...", flush=True)

# Check if run is finished sequencing
while not isRunCompleted(basePath, header["Seq_Type"]):
    print("Waiting... ({})".format(datetime.now().strftime("%H:%M:%S")), flush=True)
    time.sleep(15)#*60)

# time.sleep(15*60) # Extra wait to make sure everything is done

print("Starting pipeline launcher...", flush=True)

# Add barcodes to the respective sequencing type

allSamples[samplePosCol] = allSamples[barcodeCol]
allBarcodes = range(1,96)
if (header["Seq_Type"].lower() == "nanopore"):
    label = lambda x: f"_barcode{x}" if int(x) > 10 else f"_barcode0{x}"
    allSamples[barcodeCol] = allSamples[barcodeCol].apply(label)
    allBarcodes = [label(x) for x in allBarcodes]
elif (header["Seq_Type"].lower() == "illumina"):
    label = lambda x: f"_S{x}_"
    allSamples[barcodeCol] = allSamples[barcodeCol].apply(label)
    allBarcodes = [label(x) for x in allBarcodes]

# Split sequencing run into respective folders
print("Copying files...", flush=True)
for group in groups:
    if group.lower() == "ignore": continue

    # Check for directory
    ignore = False
    try: ignore = directories[group].lower() == "ignore"
    except KeyError: ignore = True
    if (ignore): print(f"   Ignoring file copy for {group}"); continue

    # Move files    
    outDir = directories[group]
    print("   Moving " + group + " to " + outDir, flush=True)
    includeSamples = allSamples.loc[allSamples['Sample_Group'] == group][barcodeCol].values.tolist()
    includeSamples = st.sortDigitSuffix(list(includeSamples))
    print("      Extracting samples: " + ", ".join([x.strip("_") for x in includeSamples]), flush=True)
    excludeSamples = list(set(allBarcodes) - set(includeSamples)) + ["fail","skip","unclassified","Undetermined"]
    excludeSamples = [f"*{sample}*" for sample in excludeSamples]
    shutil.copytree(basePath, outDir, ignore=shutil.ignore_patterns(*excludeSamples), dirs_exist_ok=True)

# Setup pipeline
print("Configuring pipelines...", flush=True)
for group in groups:
    if group.lower() == "ignore": continue

    # Check for pipeline
    ignore = False
    try: ignore = pipelines[group].lower() == "ignore"
    except KeyError: ignore = True
    if (ignore): print(f"   Ignoring pipeline for {group}"); continue

    # Check output dir exists
    if not os.path.exists(directories[group]):
        print(f"   Output directory does not exist for {group}. Skipping.")
        continue

    print(f"   Generating SLURM for {group}...")

    # Parse controls
    samples = allSamples.loc[allSamples['Sample_Group'] == group]
    ctrls = samples.dropna(subset=['Control'])
    samples = samples[samples['Control'].isna()]

    negCtrls = ctrls.loc[ctrls['Control'] == "negative"]
    if (len(negCtrls)):
        negCtrls = ",".join(map(str,negCtrls[samplePosCol].values.tolist()))

    posCtrls = ctrls.loc[ctrls['Control'] != "negative"]
    if (len(posCtrls)):
        posCtrls["Control"] = posCtrls[samplePosCol].astype(str) +","+ posCtrls["Control"].astype(str)
        posCtrls = " ".join(posCtrls["Control"].values.tolist())

    # Create the SLURM command
    symlink = lambda dir,link: "ln -s {} {}".format(st.findFile(os.path.join("./**",dir))[0],os.path.join(directories[group],link))

    def symlink(dir,link): # Creates symlink command
        path = os.path.join(directories[group],"**",dir)
        file = st.findFile(path)
        if (len(file) < 1): 
            raise Exception("Error: No directory found for '{}'".format(path))
        elif(len(file) > 1):
            raise Exception("Error: More than 1 directory found for '{}'".format(path))
        file = os.path.relpath(file[0],directories[group])
        link = os.path.join(directories[group],link)
        link = os.path.relpath(link,directories[group])
        return "ln -s {} {}".format(file,link)

    commands = []
    if (group == "ncov" or group == "ncov-ww"):
        # Do symlinks
        if header["Seq_Type"] == "Nanopore":
            commands.append(symlink("fast5_pass","fast5"))   
            commands.append(symlink("fastq_pass","gup_out"))   
        elif header["Seq_Type"] == "Illumina":
            commands.append(symlink("Fastq","fastq")) 

        # Go to parent dir
        parentDir = os.path.dirname(directories[group].rstrip("/")) + "/"
        commands.append("\ncd {}\n".format(parentDir))

        basecall = 1 if header["Base_Call"].lower() == "yes" else 2

        baseDir = os.path.basename(directories[group].strip("/"))

        # Run pipeline
        if (group == "ncov"):
            command = "python {} -d {} -r {} -b {}".format(pipelines[group],parentDir,baseDir,basecall)

        if (group == "ncov-ww"):
            command = "python {} -d {} -r {} -b {} -f".format(pipelines[group],parentDir,baseDir,basecall)

        if len(posCtrls): command = command + " -p {}".format(posCtrls)
        if len(negCtrls): command = command + " -c {}".format(negCtrls)
        commands.append(command)

    # Generate the SLURM file
    SLURMfile = generateSLURM(SLURM = SLURM, jobName = group+"_"+header["Run_Name"], runName = header["Run_Name"], outputDir = directories[group], command = "\n".join(commands), email = email)
    subprocess.run(["sbatch",SLURMfile,"-v"])