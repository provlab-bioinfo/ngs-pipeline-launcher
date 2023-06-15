import pandas as pd, os, searchTools as st, shutil, io, time
from configparser import ConfigParser
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
samples = getSampleSheetDataFrame(sampleSheetPath, "Samples")
directories = getSampleSheetDataVars(sampleSheetPath, "Directories")
pipelines = getSampleSheetDataVars(sampleSheetPath, "Pipelines")
header = getSampleSheetDataVars(sampleSheetPath, "Header")
basePath = header["Run_Dir"]

# Check if run is finished sequencing
while not isRunCompleted(basePath, header["Seq_Type"]):
    print("waiting...")
    time.sleep(15*60)

# Split sequencing run into respective folders
for group in samples["Sample_Group"].values.tolist():
    outDir = directories[group]
    print("Moving " + group + " to " + outDir)
    excludeSamples = samples.loc[samples['Sample_Group'] != group]["Sample_ID"].values.tolist()
    excludeSamples = ",".join(excludeSamples).split(",")
    excludeSamples = ["*"+sample+"*" for sample in excludeSamples]
    shutil.copytree(basePath, outDir, ignore=shutil.ignore_patterns(*excludeSamples))

# Create SLURM commands