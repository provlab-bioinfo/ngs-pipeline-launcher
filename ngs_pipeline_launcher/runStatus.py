import glob, itertools, os, subprocess, re
from pathlib import Path

def isRunCompleted(runPath:str):
    """Checks whether a sequencing run is completed. 
    For Illumina, checks for the 'CompletedJobInfo.xml' file. 
    For Nanopore, checks for the 'final_summary_*.txt' file.
    :param path: The path to the experiment directory
    :param seqType: The type of sequencing. Either 'Illumina' or 'Nanopore'
    :return: If complete, the list of files found. If not complete, None.
    """
    if not os.path.exists(runPath):
        return False

    # Exit if something is actively acessing any files
    # checkIfCopying = subprocess.run(['lsof', '+D', runPath], 
    #                         stdout=subprocess.PIPE, 
    #                         stderr=subprocess.PIPE)
    
    # if (checkIfCopying.returncode != 1):
    #     return None

    # Search for the target files
    completionFiles = ["final_summary_*.txt","CompletedJobInfo.xml","RunCompletionStatus.xml"]

    # Using rglob
    found = [Path(runPath).rglob(f) for f in completionFiles]
    found = [str(s) for s in list(itertools.chain.from_iterable(found))]

    # Original
    # found = [glob.glob(os.path.join(runPath,"**",f), recursive = True) for f in completionFiles]
    # found = list(itertools.chain.from_iterable(found))

    if (len(found)):
        return found
    else:
        return None