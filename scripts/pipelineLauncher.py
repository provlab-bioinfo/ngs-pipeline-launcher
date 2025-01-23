import pandas as pd, os, search_tools as st, shutil, io, time, subprocess, argparse, tempfile, pathlib, glob, re, glob, re, itertools
from configparser import ConfigParser
from datetime import datetime
from runStatus import *
pd.options.mode.chained_assignment = None  # default='warn'
os.chdir(os.path.dirname(__file__))

defaultSampleSheet = "./"
SLURM = "/nfs/APL_Genomics/apps/production/ngs-pipeline-launcher/templates/SLURM_template.batch"
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

print(f"{currentTime()} | Pipeline launcher initialized")

# Import the arguments
parser = argparse.ArgumentParser(description='APL NGS Pipeline Launcher')
parser.add_argument("-r", "--run", help="Path to the run directory. Must contain the PipelineWorksheet.xlsx.", default = defaultSampleSheet)
parser.add_argument("-e", "--email", help="Notify status alerts by e-mail.", default = None)
args = parser.parse_args()
sampleSheetPath = args.run

# print(args.email)
email = None if args.email == "None" else args.email

print(f"{currentTime()} | Checking for sequencing completion file in '{sampleSheetPath}'...", flush=True)

# Check if run is finished sequencing
while not (completionFiles := isRunCompleted(sampleSheetPath)):#isRunCompleted(basePath, header["Seq_Type"]):
    print(f"{currentTime()} |    Waiting...", flush=True)
    time.sleep(15)#*60)

print(f"{currentTime()} | Found {completionFiles}.", flush=True)

# Read data from the sample sheet
print(f"{currentTime()} | Checking for pipeline worksheet...", flush=True)
with tempfile.NamedTemporaryFile() as sampleSheet:
    os.chdir(sampleSheetPath) # TODO: This is gross
    file = st.findFiles2(os.path.join(sampleSheetPath,"**","*PipelineWorksheet*"))
    if not isinstance(file, list): file = [file]
    file = [ f for f in file if "~$" not in f ] # Handle temporary file if open in 
    if (len(file) == 0):
        raise Exception(f"No pipeline worksheet found. Please check '{sampleSheetPath}'.")
    if (len(file) > 1):
        raise Exception(f"More than one pipeline worksheet identified. Found:\n{file}.")
    
    print(f"{currentTime()} | Found {file}.", flush=True)

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
for group in groups:
    if group.lower() == "ignore": continue

    # Check for directory
    outDir = directories[group]
    ignore = False
    try: ignore = outDir.lower() == "ignore"
    except KeyError: ignore = True
    if (ignore): print(f"{currentTime()} |    Ignoring file copy for {group}"); continue

    # Move files
    print(f"{currentTime()} |    Moving {group} to {outDir}", flush=True)
    includeSamples = allSamples.loc[allSamples['Sample_Group'] == group][barcodeCol].values.tolist()
    includeSamples = st.sortDigitSuffix(list(includeSamples))
    print(f"{currentTime()} |       Extracting barcodes: " + ", ".join(st.collapseNumbers(includeSamples)), flush=True)
    excludeSamples = list(set(allBarcodes) - set(includeSamples))
    excludeSamples = st.sortDigitSuffix(list(excludeSamples))
    # print("      Excluding barcodes: " + ", ".join(st.collapseNumbers(excludeSamples)), flush=True)
    excludeSamples = [f"\/{sample}|\/.*_{sample}" for sample in excludeSamples]
    excludeSamples = excludeSamples + ["fail","skip","unclassified","Undetermined","~$"]
    excludeSamples = "|".join(excludeSamples)

    pathfilter = "**"

    if (group == "PulseNet"):
        pathfilter = "**/*fastq.gz"

    for p in glob.glob(pathfilter, recursive=True, root_dir=basePath):
        if os.path.isfile(os.path.join(basePath, p)) and not re.search(excludeSamples, p):     

            p_dest = p if (group != "PulseNet") else p[p.rindex('/')+1:]
            
            os.makedirs(os.path.join(outDir, os.path.dirname(p_dest)), exist_ok=True)
            shutil.copy(os.path.join(basePath, p), os.path.join(outDir, p_dest))
    
# Setup pipeline
print(f"{currentTime()} | Configuring pipelines...", flush=True)
for group in groups:
    if group.lower() == "ignore": continue

    # Check for pipeline
    ignore = False
    try: ignore = pipelines[group].lower() == "ignore"
    except KeyError: ignore = True
    if (ignore): print(f"{currentTime()} |    Ignoring pipeline for {group}"); continue

    # Check output dir exists
    if not os.path.exists(directories[group]):
        print(f"{currentTime()} |    Output directory does not exist for {group}. Skipping.")
        continue

    print(f"{currentTime()} |    Generating SLURM for {group}...", flush=True)

    # Parse controls
    samples = allSamples.loc[allSamples['Sample_Group'] == group]
    ctrls = samples.dropna(subset=['Control'])
    negCtrls = posCtrls = pd.DataFrame()
    samples = samples[samples['Control'].isna()]

    if(len(ctrls.index)):
        negCtrls = ctrls['Control'].str.lower() == "negative"
        negCtrls = ctrls.loc[negCtrls==True]
        negCtrls = ",".join(map(str,negCtrls[samplePosCol].values.tolist()))

        posCtrls = ctrls['Control'].str.lower() != "negative"
        posCtrls = ctrls.loc[posCtrls==True]
        posCtrls["Control"] = posCtrls[samplePosCol].astype(str) +","+ posCtrls["Control"].astype(str)
        posCtrls = " ".join(posCtrls["Control"].values.tolist())

    
    # Parse extra data
    accessions = allSamples.loc[allSamples['Sample_Group'] == group]

    # Create the SLURM command
    symlink = lambda dir,link: f"ln -s {st.findFiles2(os.path.join('./**',dir))[0]} {os.path.join(directories[group],link)}"

    def symlink(dir,link): # Creates symlink command
        path = os.path.join(directories[group],"**",dir)
        file = st.findFiles2(path)
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

    elif (group == "PulseNet"):
        parentDir = os.path.dirname(directories[group].rstrip("/")) + "/"
        baseDir = os.path.basename(directories[group].strip("/"))
        commands.append("conda activate pulsenet_analysis_pipeline")               
        commands.append(f"python {pipelines[group]} -d {parentDir} -r {baseDir}")

    elif (group == "fluA"):
        commands.append("mkdir fastq")
        commands.append("for var in {{1..101}}; do rsync -avr */Alignment_1/*Fastq/$var\_*.fastq.gz fastq/; done")
        commands.append("cd fastq")
        commands.append(f"for var in *.gz; do mv $var {header['Run_Name']}_$var ; done")
        commands.append("cd ..")
        commands.append("prog_dir=/nfs/APL_Genomics/apps/production/influenza/influenza-pipeline")
        commands.append('for x in $(find -L ./fastq -name "*R1*.fastq.gz"); do name="${x/.\/fastq\//}"; bash ${prog_dir}/generate-influenza-consensus.txt --r1 $x --r2 ${x/_R1_/_R2_} --db ${prog_dir}/influenzaDB-2022-12-08/ --outdir results --prefix ${name/_R*/}; done')
        commands.append(f"cat results/*/*.tsv | perl -ne 'BEGIN{{$h=0;}}if(/^contig/){{if($h == 0){{$h = 1; print $_;}}else{{next;}}}}else{{print $_;}}' > ./results/{header['Run_Name']}.summary.tsv")
        commands.append('echo "Job finished with exit code $? at: `date`"')     

    else:
        command = pipelines[group]

    # Generate the SLURM file
    SLURMfile = generateSLURM(SLURM = SLURM, 
                              jobName = group+"_"+header["Run_Name"], 
                              runName = header["Run_Name"], 
                              outputDir = directories[group], 
                              command = "\n".join(commands), 
                              email = email)
    out = subprocess.run(["sbatch",SLURMfile,"-v"], capture_output = True, text = True)
    # print(f"{currentTime()} |       {out.stdout}", flush=True)

print(f"{currentTime()} | All files transferred and pipeline initialized\n", flush=True)