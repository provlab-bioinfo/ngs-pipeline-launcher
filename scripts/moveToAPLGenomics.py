import os, shutil, time, argparse, glob, itertools
from datetime import datetime
from shutil import copytree, ignore_patterns
from runStatus import *

parser = argparse.ArgumentParser(description='APL NGS File mover')
parser.add_argument("-p", "--path", help="Path to the sequencing output folder")
parser.add_argument("-d", "--destination", help="Path to the experiment directory on APL Genomics")
# parser.add_argument("-t", "--type", help="'Illumina' or 'Nanopore'", required=False)
args = parser.parse_args()
path = args.path
dest = args.destination
# type = args.type

print(f"{currentTime()} | Automated file transfer tool started", flush=True)
print(f"{currentTime()} | Checking for sequencing completion file in '{path}'...", flush=True)
# Check if run is finished sequencing
while not (completionFiles := isRunCompleted(path)):
    print(f"{currentTime()} |    Waiting... ", flush=True)
    time.sleep(15)#*60)

print(f"{currentTime()} | Found {completionFiles}. Starting move...", flush=True)

# Copy all but completion files
ignoreFiles = list(set(["*" + os.path.basename(file) + "*" for file in completionFiles]))
shutil.copytree(path, dest, ignore=ignore_patterns(*ignoreFiles), dirs_exist_ok=True)

# Copy completion files
for file in completionFiles:
    if os.path.isfile(file):
        shutil.copy2(file, file.replace(path,dest))

print(f"{datetime.now().strftime('%H:%M:%S')} | Transfer Completed", flush=True)