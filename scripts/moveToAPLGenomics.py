import os, shutil, time, argparse, glob
from datetime import datetime

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

    file = ["final_summary_*.txt","CompletedJobInfo.xml"]
    if seqType.lower() == "nanopore":
        file = file[0]
    elif seqType.lower() == "illumina":
        file = file[1]

    found = [glob.glob(os.path.join(path,"**",f), recursive = True) for f in file]

    if (len(found)):
        print(f"{datetime.now().strftime('%H:%M:%S')} | Found {file} at {found}. Starting move...")
        return True
    else:
        return False

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
while not isRunCompleted(path, type):
    print(f"{datetime.now().strftime('%H:%M:%S')} | Waiting...", flush=True)
    time.sleep(15)#*60)

shutil.copytree(path, dest)

print(f"{datetime.now().strftime('%H:%M:%S')} | Transfer Completed", flush=True)

