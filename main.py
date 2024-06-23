# Created by Huy Duc - Functional Genomics Facility
# University of Colorado - Anschutz Medical Campus
'''
Use this main window to run step in the NGS_Analysis_Steps file.
Separated from all other functions for easy running.
Before running the script, make sure to move .gz files to raw folder (needs to be created)
'''

from TimeCost import timeCost, timestart
import time
from NGS_Analysis_Steps_v1_2023 import *

'''working_directory point to the main folder the files are in, in side, the 'raw' folder contains raw files'''
working_directory = '/Users/huy.duc/NGS_Sequencing_Data/Brie_and_SKBR3_NGS/Espinosa/Brie_and_SKBR3_NGS_R2'


'''the barcode_file must be in the raw folder along with the raw data file (.fastq.gz files)'''
barcode_file = 'barcodes.txt'


firstbasetokeep = 26
lastbasetokeep = 46
min_q_score = 32
'''which index do you want to align with?'''
index_name = 'brie_index'

def steps_to_run():
    # step1_unzip(working_directory)
    step2_fastQC_full(working_directory)
    # step2_get_total_readcount(working_directory)
    # step3_run_fastx_barcode_spliter(working_directory, barcode_file)
    # (working_directory)
    # step3_run_fastQC_split(working_directory)
    # step4_trimmer_location(working_directory, firstbasetokeep, lastbasetokeep)
    # step5_quality_filter_tool(working_directory, min_q_score)
    # step6_bowtie2_align(working_directory, index_name)
    # step7_comparing_output_sam_and_index(working_directory, index_name)
    # step8_making_histogram(working_directory)
    # step9_combine_output_files(working_directory, index_name)
    # step10_prep_for_mageck(working_directory)
    # '''
    # Run deseq2 in R first before running step10
    # '''
    # step11_combine_output_deseq2(working_directory, index_name)
    # checking_file()
    # fa_file_creator()
    # bowtie_index_builder()
    # printing_thing
    # check_individual_guide()

def main():
    print('----------------------', timestart(), '----------------------')
    steps_to_run()

if __name__ == '__main__':
    startTime = time.time()
    main()
    print('')
    print('---------------', timeCost(startTime), '---------------')
