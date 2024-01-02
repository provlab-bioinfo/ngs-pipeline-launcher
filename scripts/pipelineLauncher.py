import pandas as pd, os, search_tools as st, shutil, io, time, subprocess, argparse, tempfile, pathlib, glob, re, glob, re
from configparser import ConfigParser
from datetime import datetime
pd.options.mode.chained_assignment = None  # default='warn'
os.chdir(os.path.dirname(__file__))

defaultSampleSheet = "./"
SLURM = "../templates/SLURM_template.batch"
barcodeCol = "Barcode"
samplePosCol = "Sample_Pos"

currentTime = lambda: f"{datetime.now().strftime('%H:%M:%S')}"

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
    data = data.replace("[EMAIL]",  "NONE" if email is None else email)
    data = data.replace("[MAIL_TYPE]", "NONE" if email is None else "ALL")      
    data = data.replace("[RUN_DIR]", outputDir)
    outFile = os.path.join(outputDir,runName+"_SLURM.batch")
    file = open(outFile, "wt+")
    file.write(data)
    file.write("\n\n"+command)
    file.close()
    return(outFile)

print(f"Pipeline launcher started at: {currentTime()}\n")

# Import the arguments
parser = argparse.ArgumentParser(description='APL NGS Pipeline Launcher')
parser.add_argument("-r", "--run", help="Path to the run directory. Must contain the PipelineWorksheet.xlsx.", default = defaultSampleSheet)
parser.add_argument("-e", "--email", help="Notify status alerts by e-mail.", default = None)
args = parser.parse_args()
sampleSheetPath = args.run
print(args.email)
email = None if args.email == "None" else args.email

# Read data from the sample sheet
with tempfile.NamedTemporaryFile() as sampleSheet:
    file = st.findFile(os.path.join(args.run,"**","*PipelineWorksheet*"))
    file = [ f for f in file if "~$" not in f ]
    if (len(file) != 1):
        raise Exception(f"More than one PipelineWorksheet identified. Found:\n{file}")
    file = file[0]
    if pathlib.Path(file).suffix == ".xlsx":
        df = pd.read_excel(file)
        sampleSheetPath = sampleSheet.name
        df.to_csv(sampleSheetPath, index=False)
    else: sampleSheetPath = file

    header = getSampleSheetDataVars(sampleSheetPath, "Header") 
    pipelines = getSampleSheetDataVars(sampleSheetPath, "Pipelines")
    directories = getSampleSheetDataVars(sampleSheetPath, "Directories")   
    directories = {group: os.path.join(dir,header["Run_Name"]) for group, dir in directories.items()}
    allSamples = getSampleSheetDataFrame(sampleSheetPath, "Samples")

basePath = header["Run_Dir"]

# Check for appropriate inputs
if header["Seq_Type"].lower() not in ["nanopore","illumina"]:
    raise Exception(f"Seq_Type must be either 'Nanopore' or 'Illumina' in the pipeline worksheet. Found: {header['Seq_Type']}")

if header["Base_Call"].lower() not in ["yes","no"]:
    raise Exception(f"Base_Call must be either 'Yes' or 'No' in the pipeline worksheet. Found: {header['Base_Call']}")

groups = sorted(set(allSamples["Sample_Group"].dropna().values))
for group in groups:

    if group == "ignore":
        continue
    
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
    print(f"   Waiting... ({currentTime()})", flush=True)
    time.sleep(15)#*60)

# time.sleep(15*60) # Extra wait to make sure everything is done

# Add barcodes to the respective sequencing type
allSamples[samplePosCol] = allSamples[barcodeCol]
allBarcodes = range(1,1000)

if (header["Seq_Type"].lower() == "nanopore"):
    label = lambda x: f"_barcode{x}" if int(x) > 10 else f"_barcode0{x}"
elif (header["Seq_Type"].lower() == "illumina"):
    label = lambda x: f"{x}_S"
    
allSamples[barcodeCol] = allSamples[barcodeCol].apply(label)
allBarcodes = [label(barcode) for barcode in allBarcodes]

# Split sequencing run into respective folders
print("\nCopying files...", flush=True)
for group in groups:
    if group.lower() == "ignore": continue

    # Check for directory
    outDir = directories[group]
    ignore = False
    try: ignore = outDir.lower() == "ignore"
    except KeyError: ignore = True
    if (ignore): print(f"   Ignoring file copy for {group}"); continue

    # Move files
    print(f"   Moving {group} to {outDir}", flush=True)
    includeSamples = allSamples.loc[allSamples['Sample_Group'] == group][barcodeCol].values.tolist()
    includeSamples = st.sortDigitSuffix(list(includeSamples))
    print("      Extracting barcodes: " + ", ".join(st.collapseNumbers(includeSamples)), flush=True)
    excludeSamples = list(set(allBarcodes) - set(includeSamples))
    excludeSamples = st.sortDigitSuffix(list(excludeSamples))
    # print("      Excluding barcodes: " + ", ".join(st.collapseNumbers(excludeSamples)), flush=True)
    excludeSamples = [f"\/{sample}|\/.*_{sample}" for sample in excludeSamples]
    excludeSamples = excludeSamples + ["fail","skip","unclassified","Undetermined","~$"]
    excludeSamples = "|".join(excludeSamples)

    for p in glob.glob('**', recursive=True, root_dir=basePath):
        if os.path.isfile(os.path.join(basePath, p)) and not re.search(excludeSamples, p):            
            os.makedirs(os.path.join(outDir, os.path.dirname(p)), exist_ok=True)
            shutil.copy(os.path.join(basePath, p), os.path.join(outDir, p))
    
# Setup pipeline
print("\nConfiguring pipelines...", flush=True)
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

    print(f"   Generating SLURM for {group}...", flush=True)

    # Parse controls
    samples = allSamples.loc[allSamples['Sample_Group'] == group]
    ctrls = samples.dropna(subset=['Control'])
    samples = samples[samples['Control'].isna()]

    negCtrls = ctrls['Control'].str.lower() == "negative"
    if (any(negCtrls)):
        negCtrls = ctrls.loc[negCtrls]
        if (len(negCtrls)):
            negCtrls = ",".join(map(str,negCtrls[samplePosCol].values.tolist()))

    posCtrls = ctrls['Control'].str.lower() != "negative"
    if (any(posCtrls)):
        posCtrls = ctrls.loc[posCtrls]
        if (len(posCtrls)):
            posCtrls["Control"] = posCtrls[samplePosCol].astype(str) +","+ posCtrls["Control"].astype(str)
            posCtrls = " ".join(posCtrls["Control"].values.tolist())

    # Create the SLURM command
    symlink = lambda dir,link: f"ln -s {st.findFile(os.path.join('./**',dir))[0]} {os.path.join(directories[group],link)}"

    def symlink(dir,link): # Creates symlink command
        path = os.path.join(directories[group],"**",dir)
        file = st.findFile(path)
        if (len(file) < 1): 
            raise Exception(f"Error: No directory found for '{path}'")
        elif(len(file) > 1):
            raise Exception(f"Error: More than 1 directory found for '{path}'")
        file = os.path.relpath(file[0],directories[group])
        link = os.path.join(directories[group],link)
        link = os.path.relpath(link,directories[group])
        return f"ln -s {file} {link}"

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
        command = f"python {pipelines[group]} -d {parentDir} -r {baseDir} -b {basecall}"
        if (group == "ncov-ww"): command = command + " -f"

        if len(posCtrls): command = command + " -p {}".format(posCtrls)
        if len(negCtrls): command = command + " -c {}".format(negCtrls)
        commands.append(command)

    # Generate the SLURM file
    SLURMfile = generateSLURM(SLURM = SLURM, 
                              jobName = group+"_"+header["Run_Name"], 
                              runName = header["Run_Name"], 
                              outputDir = directories[group], 
                              command = "\n".join(commands), 
                              email = email)
    out = subprocess.run(["sbatch",SLURMfile,"-v"], capture_output = True, text = True)
    print(f"      {out.stdout}", flush=True)

print(f"All files transferred and pipeline initialized at: {currentTime()}\n", flush=True)