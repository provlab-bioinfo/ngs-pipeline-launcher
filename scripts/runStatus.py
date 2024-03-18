import glob, itertools, os
from datetime import datetime

def currentTime():
    return f"{datetime.now().strftime('%H:%M:%S')}"

def isRunCompleted(path:str, seqType: str = None):
    """Checks whether a sequencing run is completed. 
    For Illumina, checks for the 'CompletedJobInfo.xml' file. 
    For Nanopore, checks for the 'final_summary_*.txt' file.
    :param path: The path to the experiment directory
    :param seqType: The type of sequencing. Either 'Illumina' or 'Nanopore'
    :return: If complete, the list of files found. If not complete, None.
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
        return found
    else:
        return None

# from shutil import copytree,copy2
# def copyTree2(source, destination):
    
    # def copy2_verbose(src, dst):
#       print(f'Copying {src}')
#       copy2(src,dst)

    # copytree(source, destination, copy_function=copy2_verbose)

