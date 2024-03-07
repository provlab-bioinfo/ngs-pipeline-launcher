import os, shutil, time, argparse, glob, itertools
from datetime import datetime
from shutil import copytree, ignore_patterns

def isRunCompleted(path:str, seqType: str):
    """Checks whether a sequencing run is completed. 
    For Illumina, checks for the 'CompletedJobInfo.xml' file. 
    For Nanopore, checks for the 'final_summary_*.txt' file.
    :param path: The path to the experiment directory
    :param seqType: The type of sequencing. Either 'Illumina' or 'Nanopore'
    :return: True if complete, False if not
    """
    if not os.path.exists(path):
        return False
        #raise Exception("Run directory does not exist")

    file = ["final_summary_*.txt","CompletedJobInfo.xml"]
    if (seqType):
        if seqType.lower() == "nanopore":
            file = [file[0]]
        elif seqType.lower() == "illumina":
            file = [file[1]]

    found = [glob.glob(os.path.join(path,"**",f), recursive = True) for f in file]
    found = list(itertools.chain.from_iterable(found))

    if (len(found)):
        return True, found
    else:
        return False, found

parser = argparse.ArgumentParser(description='APL NGS File mover')
parser.add_argument("-p", "--path", help="Path to the sequencing output folder")
parser.add_argument("-d", "--destination", help="Path to the experiment directory on APL Genomics")
parser.add_argument("-t", "--type", help="'Illumina' or 'Nanopore'", required=False)
args = parser.parse_args()
path = args.path
dest = args.destination
type = args.type

print(f"{datetime.now().strftime('%H:%M:%S')} | Automated file transfer tool started", flush=True)

# Check if run is finished sequencing
completionFiles = []
while not len(completionFiles):
    isComplete, completionFiles = isRunCompleted(path, type)
    if not isComplete:
        print(f"{datetime.now().strftime('%H:%M:%S')} | Waiting...", flush=True)
        time.sleep(15)#*60)
        continue
    print(f"{datetime.now().strftime('%H:%M:%S')} | Found {completionFiles}. Starting move...")

# Copy all but completion files
shutil.copytree(path, dest, ignore=ignore_patterns('*final_summary_*.txt', '*CompletedJobInfo.xml'), dirs_exist_ok=True)

# Copy completion files
for file in completionFiles:
    if os.path.isfile(file):
        shutil.copy2(file, file.replace(path,dest))

print(f"{datetime.now().strftime('%H:%M:%S')} | Transfer Completed", flush=True)