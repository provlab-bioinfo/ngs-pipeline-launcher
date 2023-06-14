#Script for the analysis of ONT or Illumina covid19 data
#To run: python pl_ncov_auto_V5_test.py -r run_name -c ntc -b [1,2] -d dir_of_script_and_run
# i.e. python pl_ncov_auto_V5.py -r  221104_S_N_223 -c 95,96 -b 1 -d /apps/data/covid19/analysis/
#-b argument: 1 to basecall/demultiplex; 2 to skip basecalling and demultiplex
#-d argument: optional, default is cwd

# Version 2 -> Jan 12, 2022 Vince provided updates: 
#1) Remove code for modifying Artic, Pangolin and Nextclade files for inport into Bionumerics
#2) Update Guppy_basecaller, guppy_barcoder and nextflow commands
# Version 3 -> Jan 26, 2022
#1) Added code to count number of samples for Illumina fastq.
#2) Create dummy consensus fasta files when missing after freebayes nextflow pipeline
# Version 4 -> Feb 2, 2022
#1) Added ncov Dehosting step to ONT pathway after demultiplexing by guppy-plex
# 2022-09-13: VL added flag for skipping basecalling and demultiplexing.
# 2022-10-25: VL added freed v2 primer scheme for nanopolish
# Version 5 -> Nov 8, 2022
# 2022-11-08: VL removed hardcoded run dir; added argument allowing specification of run directory (default is cwd)
# Version 6 -> Dec 12, 2022
# 2022-12-12: VL upgraded to nextclade v2.8.0, guppy v6.4.2, ncov-dehoster v0.2.0; software, data, and conda env now hosted on APL_Genomics
# Added code to remove new guppy output from fastq_pass and gup_out folders.
# 2023-04-19: VL upgraded nextclade to v2.13.1.
# 2023-05-19: VL added ability to analyze data for freyja.
# 2023-06-12: AL added positive control checker

import sys
import os
import subprocess
import glob
import re
import csv
from xlsxwriter.workbook import Workbook
import pandas
from datetime import date
import shutil
import argparse

# ---------------------------------- Variable Declarations --------------------------------------------

# Variables that may be commonly modified
guppy_basecaller_num_callers = "24"
nanopolish_thread_num = "24"
nanopolish_qs_num = "16"
snakemake_core_num = "8"

# Flags for NGS type
ont_ngs_flag = 0
ill_ngs_flag = 0

artic_nanopore_nanopolish_dir = "articNcovNanopore_sequenceAnalysisNanopolish_articMinIONNanopolish"
artic_ill_cons_dir = "ncovIllumina_sequenceAnalysis_callConsensusFreebayes"
artic_ill_read_mapping_dir = "ncovIllumina_sequenceAnalysis_readMapping"
artic_ill_trim_primer_seq_dir = "ncovIllumina_sequenceAnalysis_trimPrimerSequences"
dehoster_fastq_dir = "./results/dehoster/run/fastq_pass/"
dehoster_fast5_dir = "./results/dehoster/run/fast5_pass/"
dehoster_seq_sum = "./results/dehoster/run/sequencing_summary.txt"
core_dump_dir = "./fastq_pass/guppy_basecaller-core-dump-db/"
core_dump_dir_gup = "./gup_out/guppy_barcoder-core-dump-db/"
today = date.today()

# Set up argument parser to check if arguments were passed during execution callable
parser = argparse.ArgumentParser(description='Enter values for NGS run.')

parser.add_argument('-r','-RunName',type=str, help='Name of the NGS run to process.')
parser.add_argument('-c','-NtcList',type=str, help='Comma separated list of NTC sample numbers.  Eg 4,45,96')
parser.add_argument('-b','-Basecalled',type=str, help='Type 1 to basecall and demultiplex; type 2 to skip these steps')
parser.add_argument('-d','-Directory',type=str, help='Full path of the directory of the script and run folder; default is cwd')
parser.add_argument('-f','-Freyja', action='store_true', help='add argument if data is to be analyzed by Freyja')

# Structure of the positive control 'position,refID' where refID matches the column in the current_pos_control_nextclade.xlsx file
def posCtrl(s):
    try:
        si = s.split(',')
        si[0] = int(si[0])
        return tuple(si)
    except:
        raise argparse.ArgumentTypeError("Positive controls must be given divided by commas and space, dot, or semicolon e.g.: '94,21-117-036013 95,21-129-012673'")

parser.add_argument('-p','-PosCtrList', type=posCtrl, nargs='+', help="Postive controls in the form 'position,refID' seperated by spaces Eg '94,21-117-036013 95,21-129-012673'")

args = parser.parse_args()

run_name = ""

# If argument was not passed in during execution, prompt user for the value
if(args.r == None):
    run_name = input("Enter the run name: ")
else:
    run_name = args.r

neg_control_list_str = ""

if(args.c == None):
    neg_control_list_str = input("Enter negative control sample numbers (Eg. 48,95,96): ")
else:
    neg_control_list_str = args.c


pos_control_list_str = ""

if(args.p == None):
    pos_control_list_str = input("Enter positive control sample numbers and strains (Eg. '94,21-117-036013 95,21-129-012673'): ").split(' ')
    pos_control_list_str = [posCtrl(ctrl) for ctrl in pos_control_list_str]
else:
    pos_control_list_str = args.p

# Prompts user for basecalling and demultiplexing status
skip = ""

if (args.b == None):
    skip = input("Type 1 to proceed with basecalling/demultiplexing; type 2 to skip basecalling/demultiplexing: ")
else:
    skip = args.b

# Get the run directory and set path variables
if (args.d == None):
    ngs_run_dir = os.getcwd() + "/"
else:
    ngs_run_dir = args.d

# Set variables for fb pipeline
fb_dir = "/nfs/APL_Genomics/apps/production/covid/fb/ncov2019-artic-nf"

if (args.f):
    fb_primer_scheme_ver = "V4.1"
    fb_primer_pairs_tsv = "/nfs/APL_Genomics/apps/production/covid/V4_primer/artic-ncov2019/primer_schemes/nCoV-2019/V4.1/SARS-CoV-2.primer_pairs.tsv"
    fb_primer_scheme = "SARS-CoV-2"
    ncov_primer_bed_file = "/nfs/APL_Genomics/apps/production/covid/V4_primer/artic-ncov2019/primer_schemes/nCoV-2019/V4.1/SARS-CoV-2.scheme.bed"
    fb_reference = "/nfs/APL_Genomics/apps/production/covid/V4_primer/artic-ncov2019/primer_schemes/nCoV-2019/V4.1/SARS-CoV-2.reference.fasta"
    fb_pipeline_command = ("nextflow run "+ fb_dir +" --illumina --prefix "+ run_name +" --directory ./fastq/ --composite_ref /nfs/APL_Genomics/apps/production/covid/fb/resources/composite_human_viral_reference.fna "+"--viral_contig_name MN908947.3 --scheme-directory /nfs/APL_Genomics/apps/production/covid/V4_primer/artic-ncov2019/primer_schemes " + " --ref " + fb_reference +" --bed " + ncov_primer_bed_file + " --scheme " + fb_primer_scheme + "--schemeVersion " + fb_primer_scheme_ver + " --primer_pairs_tsv "+ fb_primer_pairs_tsv)
else:
    fb_primer_scheme_ver = "freed_V2_nml"
    fb_primer_pairs_tsv = "/nfs/APL_Genomics/apps/production/covid/fb/resources/freed_primer_pairs.tsv"
    fb_primer_scheme = "nCoV-2019"
    ncov_primer_bed_file = "/nfs/APL_Genomics/apps/production/covid/artic-ncov2019/primer_schemes/nCoV-2019/freed_V2_nml/nCoV-2019.bed"
    fb_pipeline_command = ("nextflow run "+ fb_dir +" --illumina --prefix "+
run_name+" --directory ./fastq/ --composite_ref /nfs/APL_Genomics/apps/production/covid/fb/resources/composite_human_viral_reference.fna "+
"--viral_contig_name MN908947.3 --schemeRepoURL /nfs/APL_Genomics/apps/production/covid/test_V2_primer/primer-schemes/ --schemeDir primer-schemes --scheme "+
fb_primer_scheme + " --schemeVersion "+ fb_primer_scheme_ver +" --primer_pairs_tsv "+ fb_primer_pairs_tsv)

run_dir = ngs_run_dir + run_name
results_dir = ngs_run_dir + run_name + "/results/"
neg_control_barcode = "barcode96"
neg_control_name = run_name + "_" + neg_control_barcode + "_" + neg_control_barcode

# Dummy negative control run files
ntc_dummy_file_dir = "/nfs/APL_Genomics/apps/production/covid/negative_control_dummy_files/"

# File variables
qc_metric_cvs_file = run_name + ".qc.csv"
ill_qc_metric_csv_file = run_name +".qc.csv"
con_nanopolish_fasta_file = run_name +  ".cat"
con_ill_fasta_file = run_name + ".cat"
pang_lineage_file = run_name + "_lineage_report.csv"
nextclade_file = run_name + "_nextclade.tsv"
ncov_ref_genome_file = "/nfs/APL_Genomics/apps/production/covid/artic-ncov2019/primer_schemes/nCoV-2019/freed_V2_nml/nCoV-2019.reference.fasta"
if (args.f):
    ncov_ref_genome_file = "/nfs/APL_Genomics/apps/production/covid/V4_primer/artic-ncov2019/primer_schemes/nCoV-2019/V4.1/SARS-CoV-2.reference.fasta"
ncov_tools_work_snakefile = "/nfs/APL_Genomics/apps/production/covid/ncov-tools_1.9/workflow/Snakefile"
ngs_ntc_res_file = results_dir + run_name + "_" + neg_control_barcode + ".qc.csv" 
ncov_tools_summary_qc_file = results_dir+"qc_reports/"+run_name+"_summary_qc.tsv"
pangolin_lineage_file = results_dir+"lineages/"+run_name+"_lineage_report.csv"
nextclade_summary_file = results_dir+artic_nanopore_nanopolish_dir+"/"+run_name+"_nextclade.tsv"
ill_nextclade_summary_file = results_dir+artic_ill_cons_dir+"/"+run_name+"_nextclade.tsv"
neg_control_report_file = results_dir+"qc_reports/"+run_name+"_negative_control_report.tsv"
pos_control_report_file = results_dir+"qc_reports/"+run_name+"_positive_control_report.tsv"
ill_mixture_report_file = results_dir+"qc_reports/"+run_name+"_mixture_report.tsv"
pangolin_version_file = "software_allver_" + today.strftime("%Y%m%d") + ".txt"
current_pos_control_nextclade = "/nfs/APL_Genomics/apps/production/covid/current-pos-control-nextclade.csv"

# Linux commands
basecall_command = ('/nfs/APL_Genomics/apps/production/covid/guppy_6.4.2/ont-guppy/bin/guppy_basecaller --require_barcodes_both_ends '+
'-c dna_r9.4.1_450bps_hac.cfg -i fast5/ -s fastq_pass/ -r -x "auto" --gpu_runners_per_device ' + guppy_basecaller_num_callers +
' --num_callers ' + guppy_basecaller_num_callers)

demultiplex_command = ('/nfs/APL_Genomics/apps/production/covid/guppy_6.4.2/ont-guppy/bin/guppy_barcoder --require_barcodes_both_ends '+
'-i ./fastq_pass/pass -s ./gup_out/ --barcode_kits "EXP-NBD196"')

dehosting_command = ("nextflow run /nfs/APL_Genomics/apps/production/covid/ncov-dehoster_0.2.0 --nanopore --minimap2 --fastq_directory "+
"./gup_out/ --fast5_directory ./fast5/ --run_name dehoster --composite_minimap2_index /nfs/APL_Genomics/apps/production/covid/resources/composite_ref.mmi -resume")

nanopolish_primer_dir = "/nfs/APL_Genomics/apps/production/covid/ncov2019-artic-nf_freed/"

nanopolish_command = ("nextflow run "+ nanopolish_primer_dir +" --nanopolish --prefix "+
run_name + " --basecalled_fastq " + dehoster_fastq_dir + " --fast5_pass " + dehoster_fast5_dir + " --threads "+
nanopolish_thread_num + " --sequencing_summary " + dehoster_seq_sum + " --cache ~/cache_nextflow -qs " +
nanopolish_qs_num + " --schemeRepoURL /nfs/APL_Genomics/apps/production/covid/test_V2_primer/primer-schemes/ --schemeDir primer-schemes --scheme "+
"nCoV-2019 --schemeVersion freed_V2_nml")

ont_nextclade_command = ("/nfs/APL_Genomics/apps/production/covid/nextclade_2.13.1/nextclade run "+
"--input-dataset /nfs/APL_Genomics/apps/production/covid/nextclade_2.13.1/data/sars-cov-2/ "+
"--output-tsv " + nextclade_file +
" -p /nfs/APL_Genomics/apps/production/covid/nextclade_2.13.1/data/sars-cov-2/provlab_WTedit_primers_25Jul22.csv" + 
" --in-order " + con_nanopolish_fasta_file)

oicr_pipeline_command = ("nextflow run /nfs/APL_Genomics/apps/production/covid/oicr/ncov2019-artic-nf --illumina --prefix "+
run_name+" --directory ./fastq/ --composite_ref /nfs/APL_Genomics/apps/production/covid/oicr/resources/composite_human_viral_reference.fna "+
"--viral_contig_name MN908947.3 --schemeRepoURL /nfs/APL_Genomics/apps/production/covid/oicr/primer-schemes/ --schemeDir primer-schemes --scheme "+
"nCoV-2019 --schemeVersion freed_V2_nml --primer_pairs_tsv /nfs/APL_Genomics/apps/production/covid/oicr/resources/freed_primer_pairs.tsv")

oicr_ncov_tools_all_data_command = ("cp -r ../ncovIllumina_sequenceAnalysis_readMapping/* "+
"../ncovIllumina_sequenceAnalysis_trimPrimerSequences/* "+
"../ncovIllumina_sequenceAnalysis_callConsensusFreebayes/* .")

# Log File
log_file_name = run_dir+"/ncov_ngs_auto.log"
log_file = open(log_file_name, 'w',buffering=1)

# Function calls
# Function to convert a .csv file to .xlsx
def convert_csv_to_xlsx(csv_file):
    xlsx_file = csv_file.replace(".csv",".xlsx")

    # Create an XlsxWriter workbook object and add a worksheet.
    workbook = Workbook(xlsx_file,{'strings_to_numbers': True})
    worksheet = workbook.add_worksheet()

    # Create a CSV file reader.
    csv_reader = csv.reader(open(csv_file, 'rt'), delimiter=str(','))

    # Read the row data from the TSV file and write it to the XLSX file.
    for row, data in enumerate(csv_reader):
       worksheet.write_row(row, 0, data)

    workbook.close()

def convert_tsv_to_xlsx(tsv_file):
    xlsx_file = tsv_file.replace(".tsv",".xlsx")

    # Create an XlsxWriter workbook object and add a worksheet.
    workbook = Workbook(xlsx_file,{'strings_to_numbers': True})
    worksheet = workbook.add_worksheet()

    # Create a TSV file reader.
    tsv_reader = csv.reader(open(tsv_file, 'rt'), delimiter=str('\t'))

    # Read the row data from the TSV file and write it to the XLSX file.
    for row, data in enumerate(tsv_reader):
       worksheet.write_row(row, 0, data)

    workbook.close()

# Function to format the combined QC file from the artic nanopore pipeline
def generate_ont_qc_file_for_bionumerics(file,run_name,num_samples: int):
    output_file_name_sorted = file.replace(".csv","")+"_bionumeric_sort.csv"
    output_file_name = file.replace(".csv","")+"_bionumeric.csv"

    # Read in .csv file and sort by column A
    csvData = pandas.read_csv(file)
    csvData.sort_values(csvData.columns[0],axis=0,inplace=True)

    # Create the new ordered bionumerics file
    csvData.to_csv(output_file_name_sorted,index=False)

    # Read in new sorted file
    input_file = open(output_file_name_sorted, 'r')
    output_file = open(output_file_name,'w')

    lines = input_file.readlines()
    count = 0
    sample_count = 1

    barcode_num_list = [1,2,3,4,5,6,7,8,9]

    for line in lines:
        count += 1

        if(sample_count in barcode_num_list):
            barcode_num = "0"+str(sample_count)
        else:
            barcode_num = str(sample_count)

        exp_sample_name = run_name+"_barcode"+barcode_num

        line_values = line.split(",")
        sample_name = line_values[0]

        # Only write the header line as first line of new file
        if(count == 1):
            output_file.write(line.replace("\n","")+",ABPHL_key\n")
        elif(sample_name != "sample_name"):
            if(sample_name == exp_sample_name):
                output_file.write(line.replace("\n","")+","+run_name+"_"+str(sample_count)+"\n")
            else:
               # Add new row
               output_file.write(exp_sample_name+",,,,,,,,"+run_name+"_"+str(sample_count)+"\n")
               sample_count += 1
            sample_count += 1

    # Add missing samples at the end
    while(sample_count <= num_samples):
        if(sample_count in barcode_num_list):
            barcode_num = "0"+str(sample_count)
        else:
            barcode_num = str(sample_count)

        exp_sample_name = run_name+"_barcode"+barcode_num

        output_file.write(exp_sample_name+",,,,,,,,"+run_name+"_"+str(sample_count)+"\n")

        sample_count += 1

    input_file.close()
    output_file.close()

    # Convert .csv to .xslx file
    convert_csv_to_xlsx(output_file_name)

    os.remove(output_file_name_sorted)
    os.remove(output_file_name)

 
# Function to format the QC file from the artic illumina pipeline
def generate_ill_qc_file_for_bionumerics(file,run_name,num_samples: int):
    input_file = open(file,'r')
    output_file_name = file.replace(".csv","")+"_bionumeric.csv"
    output_file_name_sorted = file.replace(".csv","")+"_bionumeric_sort.csv"

    output_file = open(output_file_name, 'w')

    lines = input_file.readlines()
    count = 0

    for line in lines:
        count += 1

        line_values = line.split(",")
        sample_name = line_values[0]

        order_name = sample_name.replace(run_name+"_","")
        order_name = re.sub("_S[0-9]{1,3}_L001","",order_name)

        if(count == 1):
            output_file.write(line.replace("\n","")+",order\n")
        else:
            output_file.write(line.replace("\n","")+","+order_name+"\n")

    input_file.close()
    output_file.close()

    # Open new Bionumeric file and create a sorted list based on order column
    csvData = pandas.read_csv(output_file_name)
    csvData.sort_values(csvData.columns[8],axis=0,inplace=True)

    # Create the new ordered bionumerics file
    csvData.to_csv(output_file_name,index=False)

    # Open newly sorted file and check for missing entries
    input_file = open(output_file_name,'r')
    output_file = open(output_file_name_sorted,'w')

    lines = input_file.readlines()
    count = 0
    sample_count = 1

    for line in lines:
        count += 1

        exp_sample_name = run_name+"_"+str(sample_count)+"_S"+str(sample_count)+"_L001"

        line_values = line.split(",")
        sample_name = line_values[0]

        # Only write the header line as first line of new file
        if(count == 1):
            output_file.write(line.replace("\n","")+",ABPHL_key\n")
        elif(sample_name != "sample_name"):
            if(sample_name == exp_sample_name):
                output_file.write(line.replace("\n","")+","+run_name+"_"+str(sample_count)+"\n")
            else:
               # Add new row
               output_file.write(exp_sample_name+",,,,,,,,,"+run_name+"_"+str(sample_count)+"\n")
               sample_count += 1
            sample_count += 1

    input_file.close()
    output_file.close()

    # Convert .csv to .xslx file
    convert_csv_to_xlsx(output_file_name_sorted)

    os.remove(output_file_name)
    os.remove(output_file_name_sorted)

    xlsx_file_from = output_file_name_sorted.replace(".csv",".xlsx")
    xlsx_file_to = output_file_name.replace(".csv",".xlsx")

    os.rename(xlsx_file_from,xlsx_file_to)

def generate_pangolin_file_for_bionumerics(file):
    input_file = open(file,'r')
    output_file_name = file.replace(".csv","")+"_bionumeric.csv"
    output_file = open(output_file_name, 'w')

    lines = input_file.readlines()
    count = 0

    for line in lines:
        count += 1

        line_values = line.split(",")
        taxon_name = line_values[0]

        # Replacement for ONT sequencing
        key_value = taxon_name.replace("/ARTIC/nanopolish","")
        # Replacement for Ill sequencing
        key_value = key_value.replace("Consensus_","")

        # Only write the header line as first line of new file
        if(count == 1):
            output_file.write(line.replace("\n","")+",key\n")
        else:
            output_file.write(line.replace("\n","")+","+key_value+"\n")

    input_file.close()
    output_file.close()

    # Convert .csv to .xslx file
    convert_csv_to_xlsx(output_file_name)

    os.remove(output_file_name)

def check_neg_control_file_for_failures(file):
    input_file = open(file,'r')

    lines = input_file.readlines()
    count = 0

    for line in lines:
        count += 1

        # Skip the header
        if(count > 1):
            line_values = line.split("\t")
            neg_control = line_values[0]
            qc_status = line_values[1]

            if(qc_status != "PASS"):
                print("WARNING!!! Negative Control "+neg_control+" has a non-PASS QC status\n")
                log_file.write("WARNING!!! Negative Control "+neg_control+" has a non-PASS QC status\n")

    input_file.close()

def check_pos_control_file_for_failures(file):
    input_file = pandas.read_csv(file, sep="\t")
    past_ctrl = pandas.read_csv(current_pos_control_nextclade)
    bad_ctrl = input_file.index.to_numpy()[input_file.apply(lambda x: x['Nextclade_pango'] != x['Expected_lineage'], axis=1)].tolist()
    
    if (len(bad_ctrl)):
        for ctrl in bad_ctrl:
            print("WARNING!!! Positive Control "+ str(ctrl) +" has non-matching lineage\n")
            log_file.write("WARNING!!! Positive Control "+str(ctrl)+" has non-matching lineage\n")

    compare_cols = ["substitutions","deletions","insertions","frameShifts"]

    for index, row in input_file.iterrows():       
        known = past_ctrl.loc[past_ctrl['refID'].str.contains(row['refID'], case=False)]
        if (len(known.index) != 1):
            raise Exception("Error parsing Nextclade file.")
 
        for col in compare_cols:
            sample = set(str(row[col]).split(","))
            ref = set(str((known.iloc[0])[col]).split(","))
            diff = len(sample.symmetric_difference(ref))
            
            if (diff > 0):
                print("WARNING!!! Positive Control "+ str(index) +" differs from the expected by "+str(diff)+" "+col+"\n")
                log_file.write("WARNING!!! Positive Control "+ str(index) +" differs from the expected by "+str(diff)+" "+col+"\n")

def check_mixture_file_for_samples(file):
    input_file = open(file,'r')

    lines = input_file.readlines()
    count = 0

    for line in lines:
        count += 1

        # Skip the header
        if(count > 1):
            line_values = line.split("\t")
            sample_name_a = line_values[0]
            sample_name_b = line_values[1]

            print("WARNING!!! Samples in mixture file "+sample_name_a+" and "+sample_name_b+"\n")
            log_file.write("WARNING!!! Samples in mixture file "+sample_name_a+" and "+sample_name_b+"\n")

    input_file.close()


def AddMissingConsensusFiles(sample_list):
    print("Adding missing sample consensus fasta files")
    log_file.write("Adding missing sample consensus fasta files\n\n")

    consensus_pattern = re.compile(".*consensus.fasta")
    variant_pattern = re.compile(".*variants.norm.vcf")
    mapped_prim_trim_bam_pattern = re.compile(".*mapped.primertrimmed.sorted.bam")
    mapped_bam_pattern = re.compile(".*mapped.bam")
    sorted_bam_pattern = re.compile(".*sorted.bam")

    os.chdir(ntc_dummy_file_dir)

    # Grab list of dummy neg control run files from dir
    ntc_dummy_file_list = glob.glob("dummy.*")

    os.chdir(results_dir)

    for sample in sample_list:

        consensus_file = artic_ill_cons_dir + "/" + sample + ".consensus.fasta"
       
        if not(os.path.isfile(consensus_file)):
            print("Copying over dummy ntc files for sample: " + sample + '\n')
            log_file.write("Copying over dummy ntc files for sample: " + sample + '\n')
            
            # Iterate through ntc dummy file list and copy each to appropriate dir
            for ntc_dummy_file in ntc_dummy_file_list:
    
                from_file = ntc_dummy_file_dir + ntc_dummy_file
                sub_file = re.sub("dummy",sample,ntc_dummy_file)
                to_file = ""

                # Determine where to copy the dummy file to
                if consensus_pattern.match(ntc_dummy_file):
                    to_file = artic_ill_cons_dir + "/" + sub_file
                elif variant_pattern.match(ntc_dummy_file):
                    to_file = artic_ill_cons_dir + "/" + sub_file
                elif mapped_prim_trim_bam_pattern.match(ntc_dummy_file):
                    to_file = artic_ill_trim_primer_seq_dir + "/" + sub_file
                elif mapped_bam_pattern.match(ntc_dummy_file):
                    to_file = artic_ill_trim_primer_seq_dir + "/" + sub_file
                elif sorted_bam_pattern.match(ntc_dummy_file):
                    to_file = artic_ill_read_mapping_dir + "/" + sub_file

                if(to_file != ""):
                    print("Creating file: " + to_file + '\n')
                    log_file.write("Creating file: " + to_file + '\n')
                    shutil.copy(from_file,to_file)

                    # For consensus.fasta file, update the header to display the sample name
                    if consensus_pattern.match(ntc_dummy_file):
                        # Read in file contents to a string
                        fr = open(to_file,"r")
                        data = fr.read()
                        fr.close()

                        # Replace dummy header with sample name
                        data = re.sub("DUMMY_HEADER",sample,data)

                        # Replace contents of consensus file with updated string
                        fw = open(to_file,"w")
                        fw.write(data)
                        fw.close()

# Main Code

# Check that the Run Dir exists
if not(os.path.exists(run_dir)):
	log_file.write("The run path: "+run_dir+" does not exist, exiting\n")
	sys.exit("The run path: "+run_dir+" does not exist, exiting\n")

# Identify which type of sequencing needs to be processed
# Step 5
os.chdir(run_dir)

if(os.path.exists("fast5")):
     ont_ngs_flag = 1
elif(os.path.exists("fastq")):
     ill_ngs_flag = 1
else:
     log_file.write("Unable to determine sequence run type, exiting\n")
     sys.exit("Unable to determine sequence run type, exiting\n")


if(ont_ngs_flag == 1):
    print("Starting ONT Automated Process!\n")
    log_file.write("Starting ONT Automated Process!\n\n")
    neg_control_list = neg_control_list_str.split(',')
    
    # Iterate over list of negative controls and create a new list with full barcode names
    barcode_num_list = [1,2,3,4,5,6,7,8,9]
    barcode_neg_control_list = []
    
    for i in neg_control_list:
        barcode_num = i
        if(int(i) in barcode_num_list):
            barcode_num = "0"+str(barcode_num)

        neg_control_barcode = run_name+"_barcode"+barcode_num
        barcode_neg_control_list.append(neg_control_barcode)

    log_file.write("Negative Control Samples\n")
    # Print negative controls to logfile
    for c in barcode_neg_control_list:
        log_file.write(c+"\n")
    log_file.write('\n')

    # Steps 6 - 11
    # Manually copy the fast5 directory for the ONT sequence run

    # Step 12
    if (skip == "1"):
        print("Starting Basecalling on GPU1\n")
        log_file.write("Starting Basecalling on GPU1\n\n")

        # Check to see if fastq_pass directory exists
        if not(os.path.exists("fastq_pass")):
            subprocess.check_output(basecall_command,shell=True)

        # Move guppy basecaller core dump folder from fastq_pass
        if (os.path.exists("fastq_pass/guppy_basecaller-core-dump-db")):
            shutil.move(core_dump_dir, run_dir)

        # Step 13
        print ("Starting Demultiplex on GPU1\n")
        log_file.write("Starting Demultiplex on GPU1\n\n")

        # Ensure that the fastq_pass directory exits
        if not(os.path.exists("fastq_pass")):
            log_file.write("The fastq_pass directory is missing, exiting\n")
            sys.exit("The fastq_pass directory is missing, exiting\n")

        if not(os.path.exists("gup_out")):
            subprocess.check_output(demultiplex_command,shell=True)

    # Identify number of samples in the run
    cmd = subprocess.check_output("find gup_out/barcode* -maxdepth 1 -type d -print | wc -l",shell=True)
    num_samples = int(cmd.decode("utf-8"))

    print("Number of samples in the run = "+str(num_samples)+"\n")
    log_file.write("Number of samples in the run = "+str(num_samples)+"\n\n")

    # Move guppy basecaller core dump folder from gup_out
    if (os.path.exists("gup_out/guppy_barcoder-core-dump-db")):
            shutil.move(core_dump_dir_gup, run_dir)

    print ("Activating artic-ncov2019 Environment and Running Dehoster\n")
    log_file.write("Activating artic-ncov2019 Environment and Running Dehoster\n\n")

    # Ensure the gup_out directory exists
    if not(os.path.exists("gup_out")):
        log_file.write("The gup_out directory is missing, exiting\n")
        sys.exit("The gup_out directory is missing, exiting\n")

    if not(os.path.exists("results")):
        subprocess.run(["bash", "-c", "source activate artic-ncov2019 && "+dehosting_command+" && source exit"])

    #Check that the results directory exists
    if not(os.path.exists("results")):
        log_file.write("The results directory is missing, exiting\n")
        sys.exit("The results directory is missing, exiting\n")

# Check that new sequncing summary file was generated
    if not(os.path.isfile(dehoster_seq_sum)):
        log_file.write("Sequencing summary file was not generated, exiting\n")
        sys.exit("Sequencing summary file was not generated, exiting\n")

    # step 14 & 15
    print ("Activating artic-ncov2019 Environment and Running Nanopolish Pipeline\n")
    log_file.write("Activating artic-ncov2019 Environment and Running Nanopolish Pipeline\n\n")

    consensus_dir = results_dir + artic_nanopore_nanopolish_dir

    if not(os.path.exists(consensus_dir)):
        subprocess.run(["bash", "-c", "source activate artic-ncov2019 && "+nanopolish_command+" && source exit"])

    os.chdir(results_dir)

    # Step 17
    if not(os.path.exists(artic_nanopore_nanopolish_dir)):
        log_file.write("The "+ artic_nanopore_nanopolish_dir +" is missing, exiting\n")
        sys.exit("The "+ artic_nanopore_nanopolish_dir +" is missing, exiting\n")

    os.chdir(artic_nanopore_nanopolish_dir)

    print("Concatenating all consensus genomes into " + con_nanopolish_fasta_file+"\n")
    log_file.write("Concatenating all consensus genomes into " + con_nanopolish_fasta_file+"\n\n")

    if not(os.path.isfile(con_nanopolish_fasta_file)):
        subprocess.check_output("cat *.consensus.fasta >"+ con_nanopolish_fasta_file, shell=True)

    # Step 18
    print ("Update Pangolin DBs\n")
    log_file.write("Update Pangolin DBs\n\n")

    if(not os.path.isfile("lineage_report.csv") and not os.path.isfile(pang_lineage_file)):
        subprocess.run(["bash", "-c", "source activate pangolin && pangolin --update && exit"])

        # Grab Pangolin versions and write to file
        subprocess.run(["bash", "-c", "source activate pangolin && pangolin --all-versions > "+pangolin_version_file,"&& exit"])
        subprocess.run(["bash", "-c", "/nfs/APL_Genomics/apps/production/covid/nextclade_2.13.1/nextclade -V >> "+pangolin_version_file])

        # Step 19
        print("Running Pangolin on combined FASTA file\n")
        log_file.write("Running Pangolin on combined FASTA file\n\n")

        subprocess.run(["bash","-c", "source activate pangolin && pangolin "+ con_nanopolish_fasta_file,"&& exit"])

    # Step 20
    print("Renaming Pangolin lineage report\n") 
    log_file.write("Renaming Pangolin lineage report\n\n")

    if (os.path.isfile("lineage_report.csv")):
        subprocess.check_output("mv lineage_report.csv "+ pang_lineage_file, shell=True)

    # Step 21
    print("Running Nexclade on combined FASTA file\n")
    log_file.write("Running Nexclade on combined FASTA file\n\n")

    if not(os.path.isfile(nextclade_file)):
        subprocess.check_output(ont_nextclade_command, shell=True)

    # Step 22
    print("Moving to Results Directory\n")
    log_file.write("Moving to Results Directory\n\n")

    os.chdir(results_dir)

    # Step 23 - 27
    print("Creating config.yaml file\n")
    log_file.write("Creating config.yaml file\n\n")

    if not(os.path.isfile("config.yaml")):

        config_yaml_neg_control_line = "negative_control_samples: [ "
        config_yam_neg_controls = ""

        # Determine if any negative controls need to appear in the config file 
        for i in barcode_neg_control_list:
            nanopolish_neg_control_file = artic_nanopore_nanopolish_dir+"/"+i+".consensus.fasta"

            if(os.path.isfile(nanopolish_neg_control_file)):
               if(config_yam_neg_controls == ""):
                   config_yam_neg_controls = '\"'+ i + '\"'
               else:
                   config_yam_neg_controls = config_yam_neg_controls+',\"'+ i + '\"'

        if(config_yam_neg_controls != ""):
            config_yaml_neg_control_line = config_yaml_neg_control_line + config_yam_neg_controls + " ]"
        else:
            config_yaml_neg_control_line = "#"+config_yaml_neg_control_line + " ]"

        f = open("config.yaml", "w")
        f.write("#path to the top-level directory containing the analysis results\n")
        f.write("data_root: "+artic_nanopore_nanopolish_dir+"\n\n")
        f.write('# optionally the plots can have a "run name" prefix. If this is not defined the prefix will be "default"\n')
        f.write("run_name: "+run_name+"\n\n")
        f.write("# path to the nCov reference genome\n")
        f.write("reference_genome: "+ncov_ref_genome_file+"\n\n")
        f.write('# the sequencing platform used, can be "oxford-nanopore" or "illumina"\n')
        f.write('platform: "oxford-nanopore"\n\n')
        f.write("# path to the BED file containing the primers, this should follow the format downloaded from the ARTIC primer\n")
        f.write("primer_bed: "+ncov_primer_bed_file+"\n\n")
        f.write('# list the type of amplicon BED file that will be created from the "primer_bed".  This can include:\n')
        f.write('# full -- amplicons including primers and overlaps listed in the primer BED file\n')
        f.write('# no_primers -- amplicons including overlaps but with primers removed\n')
        f.write('# unique_amplicons -- distinct amplicons regions with primers and overlapping regions removed\n')
        f.write("bed_type: unique_amplicons\n\n")
        f.write("# minimum completeness threshold for inclusion to the SNP tree plot, if no entry\n")
        f.write("# is provided the default is set to 0.75\n")
        f.write("completeness_threshold: 0.9\n\n")
        f.write("# if a list of sample IDs for negative controls is provided, a report containing the amount\n")
        f.write("# of coverage detected in the negative controls can be generated\n")
        f.write(config_yaml_neg_control_line + '\n\n')
        f.write("# set this flag to true to include lineage assignments with pangolin in the output plots\n")
        f.write("assign_lineages: true\n")
        f.write("mutation_set: spike_mutations\n")

        f.close()

        # Step 28
        print("Activating ncov-tools environment and running snakemake\n")
        log_file.write("Activating ncov-tools environment and running snakemake\n\n")

        subprocess.run(["bash", "-c", "source activate ncov-qc_1.9 && snakemake -s "+ncov_tools_work_snakefile+ " all --cores " + snakemake_core_num + " && exit"])

    # Error Checking

    # Steps B-34 - 35
    if(os.path.isfile(neg_control_report_file)):
        print("Negative control file exists: checking for non-passing controls\n")
        log_file.write("Negative control file exists: checking for non-passing controls\n\n")
        check_neg_control_file_for_failures(neg_control_report_file)
    else: 
        print("No negative control file exists\n")
        log_file.write("No negative control file exists\n\n")

elif(ill_ngs_flag == 1):
    print("Starting Illumina Automated Process!\n")
    log_file.write("Starting Illumina Automated Process!\n\n")

    neg_control_list = neg_control_list_str.split(',')

    # Iterate over list of negative controls and create a new list with full ill names
    ill_num_list = [1,2,3,4,5,6,7,8,9]
    ill_neg_control_list = []

    for i in neg_control_list:
        ill_num = i
        if(int(i) in ill_num_list):
            ill_num = "0"+str(ill_num)

        neg_control_ill = run_name+"_"+ill_num+"_S"+ill_num+"_L001"
        ill_neg_control_list.append(neg_control_ill)
        # Add variant in naming with an extra 'S'
        neg_control_ill = run_name+"_S"+ill_num+"_S"+ill_num+"_L001"
        ill_neg_control_list.append(neg_control_ill)

    log_file.write("Negative Control Samples\n")
    # Print negative controls to logfile
    for c in ill_neg_control_list:
        log_file.write(c+"\n")
    log_file.write('\n')

    os.chdir("fastq")

    # Step 40
    print("Renaming fastq .gz files with run name\n")
    log_file.write("Renaming fastq .gz files with run name\n\n")

    # Check if fastq files have already be re-named with run-name
    cmd = subprocess.check_output('find . -maxdepth 1 -type f -name "'+run_name+'*fastq.gz" | wc -l',shell=True)
    num_runname_samples = int(cmd.decode("utf-8"))

    if(num_runname_samples == 0):
        subprocess.check_output("for x in *.gz; do mv $x "+run_name+"_$x; done",shell=True)

    # Determine the number of samples from fastq dir
    print("Determining the number of samples run\n")
    log_file.write("Determining the number of samples run\n\n")

    undetermined_pattern = re.compile(".*Undetermined.*")

    sample_list = []

    # Create a list of raw samples
    r1_sample_list = glob.glob("*R1_001.fastq.gz")
    for sample in r1_sample_list:
	    # Cutoff point for the sample name
        end_point = sample.index('_R1')

        if not undetermined_pattern.match(sample):
            sample_list.append(sample[0:end_point])

    print("Number of samples in the run = "+str(len(sample_list))+"\n")
    log_file.write("Number of samples in the run = "+str(len(sample_list))+"\n\n")

    # Step 41
    os.chdir(run_dir)

    # Step 42 + 43
    print ("Activating ill-artic-ncov2019 Environment and Running Pipeline\n")
    log_file.write("Activating ill-artic-ncov2019 Environment and Running Pipeline\n\n")
    
    if not(os.path.exists("results")):
        #subprocess.run(["bash", "-c", "source activate ill-artic-ncov2019-fb && "+oicr_pipeline_command+" && exit"])
        subprocess.run(["bash", "-c", "source activate ill-artic-ncov2019-fb && "+fb_pipeline_command+" && exit"])

    # Step 44
    #Check that the results directory exists
    if not(os.path.exists("results")):
        log_file.write("The results directory is missing, exiting\n")
        sys.exit("The results directory is missing, exiting\n")

    os.chdir(results_dir)

    if not(os.path.exists(artic_ill_cons_dir)):
        log_file.write("The "+ artic_ill_cons_dir +" is missing, exiting\n")
        sys.exit("The "+ artic_ill_cons_dir +" is missing, exiting\n")

    # Create dummy consensus FASTA files for missing samples
    if not(args.f):
        AddMissingConsensusFiles(sample_list)

    os.chdir(artic_ill_cons_dir)

    print("Concatenating all consensus genomes into " + con_ill_fasta_file+"\n")
    log_file.write("Concatenating all consensus genomes into " + con_ill_fasta_file+"\n\n")

    if not(os.path.isfile(con_ill_fasta_file)):
        subprocess.check_output("cat *.consensus.fasta >"+ con_ill_fasta_file, shell=True)

    # Step 45
    print("Update Pangolin DBs\n")
    log_file.write("Update Pangolin DBs\n\n")

    if(not os.path.isfile("lineage_report.csv") and not os.path.isfile(pang_lineage_file)):
        subprocess.run(["bash", "-c", "source activate pangolin && pangolin --update && exit"])

        # Grab Pangolin versions and write to file
        subprocess.run(["bash", "-c", "source activate pangolin && pangolin --all-versions > "+pangolin_version_file,"&& exit"])
        subprocess.run(["bash", "-c", "/nfs/APL_Genomics/apps/production/covid/nextclade_2.13.1/nextclade -V >> "+pangolin_version_file])

        # Step 46
        print("Running Pangolin on combined FASTA file\n")
        log_file.write("Running Pangolin on combined FASTA file\n\n")

        subprocess.run(["bash","-c", "source activate pangolin && pangolin "+ con_ill_fasta_file,"&& exit"])

    print("Renaming Pangolin lineage report\n") 
    log_file.write("Renaming Pangolin lineage report\n\n")

    if (os.path.isfile("lineage_report.csv")):
        subprocess.check_output("mv lineage_report.csv "+ pang_lineage_file, shell=True)

    # Step 47
    print("Running Nexclade on combined FASTA file\n\n")
    log_file.write("Running Nexclade on combined FASTA file\n")

    if not(os.path.isfile(nextclade_file)):
        subprocess.check_output(ont_nextclade_command, shell=True)

    # Step 48
    print("Create all_data fir for ncov-tools analysis\n") 
    log_file.write("Create all_data fir for ncov-tools analysis\n\n")

    os.chdir(results_dir)

    if not(os.path.exists("all_data")):
        os.makedirs("all_data")
        os.chdir("all_data")

        subprocess.check_output(oicr_ncov_tools_all_data_command,shell=True)
        os.chdir(results_dir)

    # Step 49 - 53
    print("Creating config.yaml file\n")
    log_file.write("Creating config.yaml file\n\n")

    if not(os.path.isfile("config.yaml")):

        config_yaml_neg_control_line = "negative_control_samples: [ "
        config_yam_neg_controls = ""

        # Determine if any negative controls need to appear in the config file
        for i in ill_neg_control_list:
            ill_neg_control_file = artic_ill_cons_dir+"/"+i+".consensus.fasta"

            if(os.path.isfile(ill_neg_control_file)):
               if(config_yam_neg_controls == ""):
                   config_yam_neg_controls = '\"'+ i + '\"'
               else:
                   config_yam_neg_controls = config_yam_neg_controls+',\"'+ i + '\"'

        if(config_yam_neg_controls != ""):
            config_yaml_neg_control_line = config_yaml_neg_control_line + config_yam_neg_controls + " ]"
        else:
            config_yaml_neg_control_line = "#"+config_yaml_neg_control_line + " ]"

        f = open("config.yaml", "w")
        f.write("#path to the top-level directory containing the analysis results\n")
        f.write("data_root: all_data\n\n")
        f.write('# optionally the plots can have a "run name" prefix. If this is not defined the prefix will be "default"\n')
        f.write("run_name: "+run_name+"\n\n")
        f.write("# path to the nCov reference genome\n")
        f.write("reference_genome: "+ncov_ref_genome_file+"\n\n")
        f.write('# the sequencing platform used, can be "oxford-nanopore" or "illumina"\n')
        f.write('platform: "illumina"\n\n')
        f.write("# path to the BED file containing the primers, this should follow the format downloaded from the ARTIC primer\n")
        f.write("primer_bed: "+ncov_primer_bed_file+"\n\n")
        f.write('# list the type of amplicon BED file that will be created from the "primer_bed".  This can include:\n')
        f.write('# full -- amplicons including primers and overlaps listed in the primer BED file\n')
        f.write('# no_primers -- amplicons including overlaps but with primers removed\n')
        f.write('# unique_amplicons -- distinct amplicons regions with primers and overlapping regions removed\n')
        f.write("bed_type: unique_amplicons\n\n")
        f.write("# minimum completeness threshold for inclusion to the SNP tree plot, if no entry\n")
        f.write("# is provided the default is set to 0.75\n")
        f.write("completeness_threshold: 0.9\n\n")
        f.write("# if a list of sample IDs for negative controls is provided, a report containing the amount\n")
        f.write("# of coverage detected in the negative controls can be generated\n")
        f.write(config_yaml_neg_control_line + '\n\n')
        f.write("# set this flag to true to include lineage assignments with pangolin in the output plots\n")
        f.write("assign_lineages: true\n")
        f.write("mutation_set: spike_mutations\n\n")
        f.write("# the naming convention for the bam files\n")
        f.write("# this can use the variables {data_root} (as above) and {sample}\n")
        f.write("# As per the example above, this will expand to run_200430/sampleA.sorted.bam for sampleA\n")
        f.write('bam_pattern: "{data_root}/{sample}.sorted.bam"\n\n')
        f.write("# the naming convention for the consensus sequences\n")
        f.write('consensus_pattern: "{data_root}/{sample}.consensus.fasta"\n\n')
        f.write("# the naming convention for the variants file, NF illumina runs typically use\n")
        f.write('# "{data_root}/{sample}.variants.tsv and oxford nanopore runs use "{data_root}/{sample}.pass.vcf.gz"\n')
        f.write('variants_pattern: "{data_root}/{sample}.variants.norm.vcf"\n')
        if (args.f):
            f.write('primer_prefix: "SARS-CoV-2"\n')
        f.close()

        # Step 54
        print("Activating ncov-tools environment and running snakemake\n")
        log_file.write("Activating ncov-tools environment and running snakemake\n\n")

        subprocess.run(["bash", "-c", "source activate ncov-qc_1.9 && snakemake -s "+ncov_tools_work_snakefile+ " all --cores " + snakemake_core_num + " && exit"])

      
    # Steps B-34 - 35
    if(os.path.isfile(neg_control_report_file)):
        print("Negative control file exists: checking for non-passing controls\n")
        log_file.write("Negative control file exists: checking for non-passing controls\n\n")
        check_neg_control_file_for_failures(neg_control_report_file)
    else: 
        print("No negative control file exists\n")
        log_file.write("No negative control file exists\n\n")

    # Step B-36
    if(os.path.isfile(ill_mixture_report_file)):
        print("Checking mixture file for samples\n")
        log_file.write("Checking mixture file for samples\n\n")
        check_mixture_file_for_samples(ill_mixture_report_file)
    
    # Freyja steps
    if(args.f):
        print("Running Freyja analysis\n")
        log_file.write("\nRunning freyja on samples\n")
        os.chdir(results_dir)
        os.makedirs("freyja")
        freyja_dir=results_dir+"freyja"
        os.chdir(freyja_dir)
        for bamfile in glob.glob(results_dir+"all_data/*.mapped.primertrimmed.sorted.bam"):
            shutil.copy(bamfile, freyja_dir)
        #Add freyja update here
        subprocess.run(["bash", "-c", "source activate freyja-env && for x in *.bam; do freyja variants $x --variants ${x/.mapped*/_variants} --depths ${x/.mapped*/.depths} --ref /nfs/APL_Genomics/apps/production/freyja/nCoV-2019.reference.fasta; done && exit"])
        os.makedirs("demix")
        subprocess.run(["bash", "-c", "source activate freyja-env && for x in *.tsv; do freyja demix $x ${x/_variants.tsv/.depths} --output demix/${x/_variants.tsv/_demix.tsv}; done && exit"])
        subprocess.run(["bash", "-c", "source activate freyja-env && freyja aggregate ./demix/ --output ./demix/"+run_name+"_freyja_agg.tsv && freyja plot ./demix/"+run_name+"_freyja_agg.tsv --output ./demix/"+run_name+"_freyja_agg.pdf && exit"])
        print("Freyja analysis complete\n")

if(os.path.isfile(nextclade_summary_file)):
    print("Nextclade summary file exists; checking for non-passing positive controls\n")
    log_file.write("Nextclade summary file exists; checking for non-passing positive controls\n\n")

    nextclade = pandas.read_csv(nextclade_summary_file, sep="\t")
    past_ctrl = pandas.read_csv(current_pos_control_nextclade)

    with open(pos_control_report_file, 'w') as file:
        cols = ["seqName", "Nextclade_pango", "qc.overallScore", "qc.overallStatus", "substitutions", "deletions", "insertions", "frameShifts"
,"missing"]

        posctrl = pandas.DataFrame(columns = cols)
        posctrl.insert(loc=2, column='Expected_lineage',value="")
        posctrl.insert(loc=3, column='refID',value="")

        for idx, ctrl in enumerate(pos_control_list_str):
            strain = ""
            if (ont_ngs_flag): strain = "_barcode" + str(ctrl[0])
            if (ill_ngs_flag): strain = "_S" + str(ctrl[0]) + "_"

            row = nextclade.loc[nextclade['seqName'].str.contains(strain, case=False)]
            if (len(row.index) == 0):
                print("Cannot find positive control "+str(idx)+" in the Nextclade output. Check the position number\n")
                log_file.write("Cannot find positive control "+str(idx)+" in the Nextclade output. Check the position number\n\n")
                continue
            row = row[cols]
            ref_ctrl = past_ctrl.loc[past_ctrl['refID'].str.contains(ctrl[1], case=False)]
            if (len(ref_ctrl.index) == 0):
                print("Multiple control references detected for positive control "+str(idx)+". Check the reference file.\n")
                log_file.write("Multiple control references detected for positive control "+str(idx)+". Check the reference file.\n\n")
                continue
            ref_ctrl = ref_ctrl.iloc[0]
            row.insert(2, 'Expected_lineage', ref_ctrl['Nextclade_pango'])
            row.insert(3, 'refID', ref_ctrl['refID'])
            posctrl = pandas.concat([posctrl,row], ignore_index=True)
        posctrl.to_csv(pos_control_report_file,"\t")
    check_pos_control_file_for_failures(pos_control_report_file)    
else: 
    print("No Nextclade summary file exists; cannot create the positive control file\n")
    log_file.write("No positive control file was generated\n\n")

log_file.close()

os.chdir(ngs_run_dir)

