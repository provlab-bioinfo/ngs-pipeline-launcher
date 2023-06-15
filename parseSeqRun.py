import pandas as pd, os, searchTools as st, shutil, io, time
from configparser import ConfigParser
pd.options.mode.chained_assignment = None  # default='warn'
os.chdir(os.path.dirname(__file__))

def getSampleSheetDataVars(path:str, section:str):
    cfg = ConfigParser(allow_no_value=True)
    cfg.optionxform = str
    cfg.read(path)
    dict = {k:v for k, *v in map(lambda x: str.split(x,sep=","), cfg[section])}
    dict = {k:v[0] for k, v in dict.items() if v} # Removes blank keys
    return (dict)

def getSampleSheetDataFrame(path:str, section:str):
    cfg = ConfigParser(allow_no_value=True)
    cfg.optionxform = str
    cfg.read(path)
    buf = io.StringIO()
    buf.writelines('\n'.join(row.rstrip(',') for row in cfg[section]))
    buf.seek(0)
    return (pd.read_csv(buf))

def readSampleSheet(path:str):
    sampleSheet = pd.read_csv(path)
    # sampleSheet['Sample_ID'] = sampleSheet['Sample_ID'].transform(lambda x: '"'+x+'"')
    sampleSheet = sampleSheet.groupby(['Sample_Group'], as_index=False).agg({'Sample_ID': ','.join})
    return(sampleSheet)

def isRunCompleted(path:str, seqType: str):
    file = ""
    if seqType == "Nanopore":
        file = "final_summary_*.txt"
    elif seqType == "Illumina":
        file = "CompletedJobInfo.xml"
    return False if not len(st.findFile(os.path.join(path,"**",file))) else True

sampleSheetPath = "SampleSheet.csv"
allSamples = getSampleSheetDataFrame(sampleSheetPath, "Samples")
directories = getSampleSheetDataVars(sampleSheetPath, "Directories")
pipelines = getSampleSheetDataVars(sampleSheetPath, "Pipelines")
header = getSampleSheetDataVars(sampleSheetPath, "Header")
basePath = header["Run_Dir"]

# Check if run is finished sequencing
while not isRunCompleted(basePath, header["Seq_Type"]):
    print("waiting...")
    time.sleep(15*60)

# Split sequencing run into respective folders
for group in set(allSamples["Sample_Group"].values):
    outDir = directories[group]
    print("Moving " + group + " to " + outDir)
    excludeSamples = allSamples.loc[allSamples['Sample_Group'] != group]["Sample_ID"].values.tolist()
    excludeSamples = ",".join(excludeSamples).split(",")
    excludeSamples = ["*"+sample+"*" for sample in excludeSamples]
    shutil.copytree(basePath, outDir, ignore=shutil.ignore_patterns(*excludeSamples))

# Create SLURM commands
for group in set(allSamples["Sample_Group"].values):
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
    if (group == "ncov"):
        command = [pipelines[group],"-d",directories[group],"-r",header["Run_Name"],"-b","2"]
        if len(posCtrls): command = command + ["-p",posCtrls]
        if len(negCtrls): command = command + ["-c",negCtrls]
        print(" ".join(command))

    elif (group == "ncov-ww"):
        command = [pipelines[group],"-d",directories[group],"-r",header["Run_Name"],"-b","2","-f"]
        if len(posCtrls): command = command + ["-p",posCtrls]
        if len(negCtrls): command = command + ["-c",negCtrls]
        print(" ".join(command))

