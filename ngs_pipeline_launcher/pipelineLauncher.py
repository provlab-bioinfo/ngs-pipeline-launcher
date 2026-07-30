import pandas as pd, os, search_tools as st, shutil, io, time, subprocess, argparse, tempfile, pathlib, re, glob
from configparser import ConfigParser
from itertools import chain
from pathlib import Path
import openpyxl as xl
from runStatus import *
os.chdir(os.path.dirname(__file__))

defaultSampleSheet = "./"
SLURM = "/nfs/APL_Genomics/apps/production/ngs-pipeline-launcher/templates/SLURM_template.batch"
barcodeCol = "Barcode" # The column in the Pipeline Worksheet [Samples] section that contains the barcode
samplePosCol = "Sample_Pos" # Generated column in the Pipeline Worksheet [Samples] section that contains the sample name (e.g., "27_S" for illumina)
symlinkFQ = False # Should fastq's be symlinked?
pathfilter = ["**/*fastq.gz","**/*fq.gz","**/*fast5","**/report_*.json","**/*PipelineWorksheet*"]  # Which files should be transferred?

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

def subsetWorksheet(group:str, path:str, outPath:str, maxCols:int = 26):
    """Subsets a worksheet to only include a specific Sample_Group in the [SAMPLES] section.
    Due to limitations with openpyxl, formatting cannot be removed from the entire row without a significant amount of computation,
    so formatting will only be changed for columns 1 to 'maxCols'.
    :param group: The group to look for in the column 'Sample_Group'
    :param path: The input path of the pipeline worksheet
    :param out_path: The output path of the subsetted pipeline worksheet
    :param maxCols: The maximum number of columns to display/format
    :return: A DataFrame representing the sections
    """ 
    wb = xl.load_workbook(path)
    ws = wb.active

    # Sets the max visible columns
    last_col = maxCols 
    for col_idx in range(last_col+1, 16385):
        col_letter = xl.utils.get_column_letter(col_idx)
        # if (ws.column_dimensions[col_letter].hidden): break
        ws.column_dimensions[col_letter].hidden = True

    # Remove samples from the [Samples] section
    rows = list(ws.iter_rows(min_row=1, max_row=ws.max_row))

    for row in reversed(rows): 
        cell = row[2] # col idx 3 is Sample_Group, TODO: Search for this instead of hardcoding
        if cell.value == "Sample_Group":
            break
        if cell.value != group:
            ws.delete_rows(idx = cell.row)
            # Clear styles from the last row in the sheet, as the data will have shifted upwards
            for row in ws.iter_cols(min_row = ws.max_row+1, min_col = 1, max_col = last_col+1, max_row = ws.max_row+1):
                for cell in row:
                    cell.style = "Normal"

    # Hide directories in the [Directories section]
    hideRowsExceptGroup(ws = ws, col = 1, section = "[Directories]", group = group)
    hideRowsExceptGroup(ws = ws, col = 1, section = "[Pipelines]",   group = group)

    wb.save(outPath)

def hideRowsExceptGroup(ws, col, section, group):
    """Hides rows in the worksheet that do not match the group
    :param ws: The worksheet to check
    :param path: The column to check, as an integer. E.g., 1 for 'A'
    :param section: The section header. Typically '[Directories]' or '[Pipelines]'.
    :param group: The group to check for (case sensitive).
    :return: A DataFrame representing the sections
    """     
    sectionIdx = 0

    # Find the start of the section
    for cell in ws[xl.utils.get_column_letter(1)]:
        if cell.value == section:
            sectionIdx = cell.row
            break
    
    for row in range(sectionIdx+1,ws.max_row+1):
        val = ws.cell(row,1).value
        if val is None: # Break if the first column is blank
            break
        val = ws.cell(row,col).value
        if val != group: # Check if the target column contains the group
            ws.row_dimensions[row].hidden = True

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

def printLog (message: str) :
    """Prints a message formatted as "[Current time] | [Message]
    :param message: The message to print
    """
    print (f"{currentTime()} | {message}", flush=True)

def runLauncher(sampleSheetPath: str, email: str = None, if_exists = "error"):
    """Generates a SLURM command file based on a template
    :param sampleSheetPath: The path to the directory containing the sample sheet. The file must contain the sub-string 'PipelineWorksheet'.
    :param email: If desired, SLURM will send an e-mail on job status. Default = None
    :param if_exists: What to do if directory already exists? Options: 'error' out , 'ignore' the group, or 'delete' the existing directory. Default = 'error'.
    """
    printLog(f"Pipeline launcher initialized for '{sampleSheetPath}'...")

    # Read data from the sample sheet
    printLog(f"Checking for pipeline worksheet...")
    with tempfile.NamedTemporaryFile() as sampleSheet:

        # Check for either specific file or directory to search
        if os.path.isfile(sampleSheetPath): # If file
            file = sampleSheetPath  
        elif os.path.isdir(sampleSheetPath): # If directory, then search for worksheet
            file = st.findFiles2(os.path.join(sampleSheetPath,"**","*PipelineWorksheet*"))
            if not isinstance(file, list): file = [file]
            file = [ f for f in file if "~$" not in f ] # Exclude temp files (e.g., if open in Excel)
            if (len(file) == 0): # If no files found
                raise Exception(f"No pipeline worksheet found. Filename must include 'PipelineWorksheet'. Please check '{sampleSheetPath}'.")
            if (len(file) > 1): # If more than 1 file is found
                raise Exception(f"More than one pipeline worksheet identified. Only 1 filename can contain 'PipelineWorksheet'. Found:\n{file}.")
            file = file[0]
        else:
            raise Exception(f"Pipeline worksheet not found. Filename must contain 'PipelineWorksheet'. Please check '{sampleSheetPath}'.")
        
        printLog(f"   Found '{file}'")
        
        # Convert to CSV
        if pathlib.Path(file).suffix.lower() == ".xlsx":
            df = pd.read_excel(file)
            sampleSheetPath = sampleSheet.name # Export df to the tempfile
            df.to_csv(sampleSheetPath, index=False)
        elif pathlib.Path(file).suffix.lower() == ".csv": 
            sampleSheetPath = file
        else:
            raise Exception(f"Pipeline worksheet must have the file type of '.xlsx' or '.csv'. Please check '{file}'.")

        # Get the variables for the run
        header = getSampleSheetDataVars(sampleSheetPath, "Header") 
        runName = header["Run_Name"].strip()
        runDir = header["Run_Dir"].strip()
            
        pipelines = getSampleSheetDataVars(sampleSheetPath, "Pipelines")
        directories = getSampleSheetDataVars(sampleSheetPath, "Directories")   
        directories = {group: os.path.join(dir,runName) for group, dir in directories.items()} # Adds path name to end of directory path
        allSamples = getSampleSheetDataFrame(sampleSheetPath, "Samples")

    # Check if run is finished sequencing
    if not os.path.isdir(runDir): # Check if run exists
        raise Exception(f"Run directory does not exist at '{header['Run_Dir']}'.")

    printLog(f"Checking for sequencing completion file...")
    sleep_time = 60
    while not (completionFiles := isRunCompleted(runDir)):#isRunCompleted(basePath, header["Seq_Type"]):
        printLog(f"   Waiting...")
        time.sleep(sleep_time)
        sleep_time = min(3600, sleep_time*2) # Use an increasing wait timer for sleeping, to max of 1 hour

    if any(".xml" in i for i in completionFiles): # "CompletedJobInfo.xml","RunCompletionStatus.xml"
        platform = "illumina"
    elif any(".txt" in i for i in completionFiles): # "final_summary_*.txt"
        platform = "nanopore"
    else:
        raise Exception(f"Cannot detect which run type (Illumina or ONT). Should only have either one of 'CompletedJobInfo.xml' and 'RunCompletionStatus.xml'for Illumina, or 'final_summary_*' for ONT in '{completionFiles}'")

    printLog(f"   Found '{completionFiles}'")

    # Check for appropriate inputs
    groups = sorted(set(allSamples["Sample_Group"].dropna().values))
    for group in groups[:]:

        # Remove group if name is 'ignore'
        if group == "ignore":
            groups.remove(group)
            continue

        # Remove group if target directory is 'ignore'
        ignore = False
        try: ignore = directories[group].lower().strip() == "ignore"
        except KeyError: ignore = True
        if (ignore): 
            printLog(f"   Ignoring file copy for {group}"); 
            groups.remove(group)
            continue
        
        # Check if the pipeline exists
        pipeline, *args = pipelines[group].split(" ")
        if not os.path.exists(pipeline):
            if (pipeline != "ignore"):
                raise Exception(f"Pipeline for '{group}' does not exist at '{pipeline}'")

        # Check if the directories exist
        if os.path.exists(directories[group]):
            if (directories[group] != "ignore"):
                if (len(directories[group]) != 0): # If directory already exists, check if_exists argument
                    if if_exists == "error":
                        raise Exception(f"Directory for '{group}' at '{directories[group]}' already exists and is not empty. Please choose empty or non-existing directory.")
                    if if_exists == "delete":
                        printLog(f"Removing dir: {directories[group]}.")
                        shutil.rmtree(directories[group], ignore_errors=True)
                    if if_exists == "ignore":
                        groups.remove(group) # Remove the group from further processing
                    
    time.sleep(5) # Extra little wait to make sure everything is done

    # Add barcodes to the respective sequencing type
    allSamples[samplePosCol] = allSamples[barcodeCol]
    allBarcodes = range(1,500) # Max number is arbitrary. Only needs to be higher than the max barcode value.

    if platform == "illumina":
        label = lambda x: f"{x}_S"
    elif platform == "nanopore":
        label = lambda x: f"barcode{x}" if int(x) >= 10 else f"barcode0{x}"
        
    allSamples[barcodeCol] = allSamples[barcodeCol].apply(label)
    allBarcodes = [label(barcode) for barcode in allBarcodes]

    # Split sequencing run into respective folders
    printLog(f"Locating for files to move...")

    for group in groups:

        # Check for directory
        outDir = directories[group]

        # Get barcodes to include
        printLog(f"   Moving {group} to '{outDir}'")
        includeSamples = allSamples.loc[allSamples['Sample_Group'] == group][barcodeCol].values.tolist()
        includeSamples = st.sortDigitSuffix(list(includeSamples))
        printLog(f"      Extracting barcodes: " + ", ".join(st.collapseNumbers(includeSamples)))

        # Get barcodes to exclude
        excludeSamples = list(set(allBarcodes) - set(includeSamples))
        excludeSamples = st.sortDigitSuffix(list(excludeSamples))
        # print("      Excluding barcodes: " + ", ".join(st.collapseNumbers(excludeSamples)))
        excludeSamples = [f"\/{sample}|\/.*_{sample}|\/.*-{sample}" for sample in excludeSamples]
        excludeSamples = excludeSamples + ["fail","skip","unclassified","Undetermined","~\$","pod5"]
        excludeSamples = "|".join(excludeSamples)

        # Move files
        found_files = list(chain.from_iterable(Path(runDir).rglob(pattern) for pattern in pathfilter))
        found_files = [str(file.relative_to(runDir)) for file in found_files if not re.search(excludeSamples, str(file))] # Remove exluded samples
        
        if not any(file.endswith(("fastq.gz","fq.gz")) for file in found_files):
            printLog(f"No fastq's found for {group}. Ignoring the pipeline.")
            pipelines[group] = "ignore"
            continue

        for filePath in found_files:
            os.makedirs(os.path.join(outDir, os.path.dirname(filePath)), exist_ok=True)
            src = os.path.join(runDir, filePath)
            dst = os.path.join(outDir, filePath)
            if symlinkFQ and pathlib.Path(filePath).suffix.lower() in [".fastq.gz", ".fq.gz"]: # Symlink or not
                os.symlink(src, dst)
            elif "PipelineWorksheet" in filePath:
                subsetWorksheet(group, src, dst)
            else:
                shutil.copy(src, dst)

        printLog(f"      Copied files: {len(found_files)}")
                
    # Setup pipeline
    printLog(f"Configuring pipelines...")
    for group in groups:

        # Check for pipeline
        ignore = False
        try: ignore = pipelines[group].lower() == "ignore"
        except KeyError: ignore = True
        if (ignore): printLog(f"   Ignoring pipeline for {group}"); continue

        # Check output dir exists
        if not os.path.exists(directories[group]):
            printLog(f"   Output directory does not exist for {group}. Skipping.")
            continue

        printLog(f"   Generating SLURM for {group}...")

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

        # Create the SLURM file
        commands = [f"cd {directories[group]}"]

        if (group != ""):

            # Check for pipeline extension to determine interpreter
            pipeline, *args = pipelines[group].split(" ")            
            if pipeline.lower().endswith((".py")):
                 type = "python"
            elif pipeline.lower().endswith((".sh")):
                 type = "bash"
            else: 
                raise Exception(f"Error: Unsure how to start generic pipeline '{pipelines[group]}'. This currently only supports '.py' and '.sh' scripts.")

            # Create commands
            command = f"{type} {pipelines[group]} -r {directories[group]}"
            if len(posCtrls): command = f"{command} -p '{posCtrls}'"
            if len(negCtrls): command = f"{command} -c {negCtrls}"
            commands.append(command)

        else:
            printLog(f"   No pipeline found for {group}. Skipping.")
            continue

        # Generate the SLURM file
        SLURMfile = generateSLURM(SLURM = SLURM, 
                                jobName = group+"_"+runName, 
                                runName = runName, 
                                outputDir = directories[group], 
                                command = "\n".join(commands), 
                                email = email)
        out = subprocess.run(["sbatch",SLURMfile,"-v"], capture_output = True, text = True) # Launch the SLURM file
        printLog(f"      {out.stdout}")

    printLog(f"All files transferred and pipeline initialized\n")

# Import the arguments
def lower_and_strip(value):
    return str(value).strip().lower()

parser = argparse.ArgumentParser(description='APL NGS Pipeline Launcher')
parser.add_argument("-r", "--run", help="Path to the run directory. Must contain the PipelineWorksheet.xlsx.", default = defaultSampleSheet)
parser.add_argument("-e", "--email", help="Notify status alerts by e-mail.", default = None)
parser.add_argument("-x", "--if_exists", help="What to do if directory already exists? Options: 'error' out , 'ignore' the group, or 'delete' the existing directory. Default: 'error'.", default = 'error', type = lower_and_strip)
args = parser.parse_args()

if (args.if_exists not in ["error",'ignore','delete']):
    raise Exception(f"Argument --if_exists must be either 'error', 'ignore', or 'delete'.")

runLauncher(args.run, None if args.email == "None" else args.email, args.if_exists)