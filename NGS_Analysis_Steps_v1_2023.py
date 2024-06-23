# Created by Huy Duc - Functional Genomics Facility
# University of Colorado - Anschutz Medical Campus
'''
All Functions Go Here!!!
'''


import os
import sys
from bioinfokit.analys import fastq
from collections import Counter, defaultdict
import shutil
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


'''Required software:'''
#bowtie2
bowtie2_folder = '/Users/huy.duc/binary/bowtie2-2.4.1-macos-x86_64'
#fastQC
fastqprogram = '/Users/huy.duc/binary/FastQC/fastqc'
#fastx_toolkit
fastx_toolkit = '/Users/huy.duc/binary/fastx_toolkit_0.0.13'
# index_collection is where the index (created by bowtie2 are located
index_collection = '/Users/huy.duc/binary/bowtie2-2.4.1-macos-x86_64/Index_collections'
# index_file: the location of fa files of the index
index_folder = '/Users/huy.duc/binary/bowtie2-2.4.1-macos-x86_64/Index_files'


def get_reverse_compliment(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    letters = list(s[::-1])
    letters = [basecomp.get(base, 'N') for base in letters]
    return ''.join(letters)


def bowtie_index_builder():
    '''Function to use bowtie as command line to generate the index for CSPC oligos'''
    bowtie2_build = bowtie2_folder + '/bowtie2-build'
    fa_file = ''  # input fa file here
    index_file = fa_file.split('/')[-1].split('.fa')[0]
    command_line = bowtie2_build + ' ' + fa_file + ' ' + bowtie2_folder + index_file
    os.system(command_line)


def step1_unzip(working_directory):
    file_location = working_directory + '/raw'
    for file in os.listdir(file_location):
        if file.endswith('.fastq.gz'): # only unzip the forward file R1, the reverse R2 is not needed
            print('Unzipping file:',file)
            file_path = file_location+'/'+file
            unzipped_file_name = file_path.split('.gz')[0]
            with gzip.open(file_path, 'rb') as f_in:
                with open(unzipped_file_name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)


def step2_fastQC_full(working_directory):
    foldernametocreate = 'Fast_QC_full'
    MYDIR = working_directory + '/raw/' + foldernametocreate
    CHECKFOLDER = os.path.isdir(MYDIR)
    if not CHECKFOLDER:
        os.makedirs(MYDIR)
        print('Created Folder', foldernametocreate)
    else:
        print('Folder', foldernametocreate, 'already exists!')


    argument = ' --nogroup '  #this will not group the nucleotides together
    for file in os.listdir(working_directory + '/raw'):
        if file.endswith('.fastq.gz'):
            print(file)
            commandline = fastqprogram + argument + working_directory + '/raw/' +file +  ' -o ' + MYDIR
            print(commandline)
            os.system(commandline)


def step2_get_total_readcount(working_directory):
    for file in os.listdir(working_directory + '/raw'):
        if file.endswith('.fastq'):
            print(file)
            total_read_count = 0 # this should be 4 times the number of reads
            read_counter = 0
            with open(working_directory+'/raw/'+file) as reader:
                for line in reader:
                    line = line.strip('\n')
                    if '+' in line or '@' in line or 'I' in line:
                        total_read_count += 1
                        continue
                    else:
                        total_read_count += 1
                        read_counter += 1
            print('Number of line in the file: {:,}'.format(total_read_count))
            print('Number of read printed: {:,}'.format(read_counter))


def step3_get_read_count(working_directory):
    '''Function to get the read count for each NGS file received'''
    for file in os.listdir(working_directory + '/Barcode_Split'):
        if file.endswith('.fastq'):
            print(file)
            read_counter = 0
            with open(working_directory + '/Barcode_Split/' + file, 'r') as reader:
                for line in reader:
                    line = line.strip('\n')
                    if '+' in line or '@' in line or 'F' in line:
                        continue
                    else:
                        read_counter += 1
            print('Number of line printed: {:,}'.format(read_counter))


def step3_run_fastx_barcode_spliter(working_directory, barcode_file):
    foldernametocreate = 'Barcode_Split'
    MYDIR = working_directory + '/' + foldernametocreate
    CHECKFOLDER = os.path.isdir(MYDIR)
    if not CHECKFOLDER:
        os.makedirs(MYDIR)
        print('Created Folder', foldernametocreate)
    else:
        print('Folder', foldernametocreate, 'already exists!')

    if barcode_file == 'none':
        pass
    else:
        barcode_file_path = working_directory + '/raw/'+barcode_file
        fastx_splitter = fastx_toolkit + '/fastx_barcode_splitter.pl'
        mismatchnumber = 1
        for file in os.listdir(working_directory + '/raw'):
            if file.endswith('.fastq') and 'R1' in file:
                fastq_file_path = working_directory + '/raw/'+file
                prefix = fastq_file_path.split('.fastq')[0] + '_'
                commandtorun = 'cat ' + fastq_file_path + ' | ' + fastx_splitter + ' --bcfile ' + barcode_file_path + ' --bol --mismatches ' + str(mismatchnumber) + ' --prefix ' + prefix + ' --suffix ".fastq"'
                print(commandtorun)
                os.system(commandtorun)

    for file in os.listdir(working_directory+'/raw/'):
        if file.endswith('.fastq'):
            shutil.move(working_directory+'/raw/'+file, MYDIR + '/' + file)


def step3_run_fastQC_split(working_directory):
    '''Function to run fastQC on split reads based on barcodes'''
    foldernametocreate = 'Fast_QC'
    MYDIR = working_directory + '/Barcode_Split/' + foldernametocreate
    CHECKFOLDER = os.path.isdir(MYDIR)
    if not CHECKFOLDER:
        os.makedirs(MYDIR)
        print('Created Folder', foldernametocreate)
    else:
        print('Folder', foldernametocreate, 'already exists!')

    argument = ' --nogroup '
    for file in os.listdir(working_directory + '/Barcode_Split'):
        if file.endswith('.fastq'):
            print(file)
            commandline = fastqprogram + argument + working_directory + '/Barcode_split/' + file + ' -o ' + MYDIR
            print(commandline)
            os.system(commandline)


def step4_trimmer_location(working_directory, firstbasetokeep, lastbasetokeep):
    '''Function to trim each read based on location, written by Huy Duc'''
    foldernametocreate = 'Trimmed_Split_Files_Location'
    MYDIR = working_directory + '/' + foldernametocreate
    CHECKFOLDER = os.path.isdir(MYDIR)
    if not CHECKFOLDER:
        os.makedirs(MYDIR)
        print('Created Folder', foldernametocreate)
    else:
        print('Folder', foldernametocreate, 'already exists!')
    fastq_file_input_directory = os.listdir(working_directory + '/Barcode_Split/')
    for fastq_file_input in fastq_file_input_directory:
        if fastq_file_input.endswith('.fastq'):
            fastq_file_output = working_directory + '/'+foldernametocreate + '/' +fastq_file_input.split('.fastq')[0] + '_trimmed_loc_'+str(firstbasetokeep)+'_' + str(lastbasetokeep) +'.fastq'
            print('Trimming file: ', fastq_file_input, 'from', firstbasetokeep, '-', lastbasetokeep)
            fastq_iter = fastq.fastq_reader(working_directory + '/Barcode_Split/'+fastq_file_input)
            with open(fastq_file_output, 'w+') as output_writer:
                for element in fastq_iter:
                    output_writer.write(element[0] + '\n')
                    output_writer.write(element[1][firstbasetokeep:lastbasetokeep] + '\n')
                    output_writer.write(element[2] + '\n')
                    output_writer.write(element[3][firstbasetokeep:lastbasetokeep] + '\n')


def step5_quality_filter_tool(working_directory, min_q_score):
    '''Function to filter read based on Q quality score, default is 32'''
    '''This function will also remove any read less than 15 bp'''
    foldernametocreate = 'Trimmed_Split_Files_' + 'location_Filtered'
    MYDIR = working_directory + '/' + foldernametocreate
    CHECKFOLDER = os.path.isdir(MYDIR)
    if not CHECKFOLDER:
        os.makedirs(MYDIR)
        print('Created Folder', foldernametocreate)
    else:
        print('Folder', foldernametocreate, 'already exists!')

    for file in os.listdir(working_directory + '/Trimmed_Split_Files_location/'):
        if file.endswith('.fastq'):
            input_path = working_directory + '/Trimmed_Split_Files_location/' + file
            print('Filtering file: ', file)
            read_fastq = fastq.fastq_reader(input_path)
            original_counter = 0
            filtered_counter = 0
            output = file.split('.fastq')[0] + '_filtered.fastq'
            output_path = working_directory + '/Trimmed_Split_Files_location_Filtered/' + output
            with open(output_path, 'w+') as outputwriter:
                for read in read_fastq:
                    original_counter += 1
                    if len(read[3]) >= 15:
                        sum_q_score = 0
                        for char in read[3]:
                            individual_score = ord(char)-33
                            sum_q_score = sum_q_score + individual_score
                        average_q_score = sum_q_score/len(read[3])
                        if average_q_score >= min_q_score:
                            filtered_counter += 1
                            outputwriter.write(read[0] + '\n')
                            outputwriter.write(read[1] + '\n')
                            outputwriter.write(read[2] + '\n')
                            outputwriter.write(read[3] + '\n')
                        else:
                            continue
                    else:
                        continue
                print('             Number of read in original file: {:,}'.format(original_counter))
                print('             Number of read in filtered file: {:,}'.format(filtered_counter))


def step6_bowtie2_align(working_directory, index_name):
    '''Function to use bowtie2 to aling fastq file to index'''

    bowtie2 = bowtie2_folder + '/bowtie2'

    trimmed_file_directory = working_directory + '/Trimmed_Split_Files_location_Filtered/'
    folder_to_create = 'Bowtie2_Aligned_Filtered_files'
    MYDIR = working_directory + '/' + folder_to_create
    CHECK_FOLDER = os.path.isdir(MYDIR)
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        print('Created folder:', folder_to_create)
    else:
        print(folder_to_create, 'folder already exists.')
    aligned_folder = working_directory + '/' + folder_to_create
    ### This block of code opens a standard output to write the output in the screen to a text file
    stdOutputFile = aligned_folder + '/stdoutput.txt'
    stdoutOrigin = sys.stdout
    sys.stdout = open(stdOutputFile, 'w')

    for file in os.listdir(trimmed_file_directory):
        if file.endswith('.fastq'):
            print('Running bowtie2 on file:', file)
            print('Index to run: ', index_name)
            file_path = trimmed_file_directory + '/' + file
            bam_output = aligned_folder + '/' + file.split('.fastq')[0] + '_aligned_with_' + index_name + '.bam'
            index_location = index_collection + '/' + index_name + '/' + index_name
            bowtie2_command = bowtie2 + ' -x ' + index_location + ' ' + file_path + ' -S ' + bam_output
            os.system(bowtie2_command)


def step7_comparing_output_sam_and_index(working_directory, index_name):
    '''Function to compare the matched output from sam file and index file'''
    for file in os.listdir(working_directory + '/Bowtie2_Aligned_Filtered_files'):
        if file.endswith('.bam'):
            ### a file that contains all of the guide not in found in NGS, just the guide ID from brunello index file
            bamfilename = file.split('/')[-1]
            print(f'Comparing {bamfilename} with {index_name}')
            index_file_path = index_folder + index_name + '.fa'
            totalCounterIndex = 0
            indexGuideCounter = 0
            guideIDList = list() ### list to hold all of guide name (start with '>')
            with open(index_file_path, 'r') as indexReader:
                for indexLine in indexReader:
                    indexLine = indexLine.strip('\n')
                    totalCounterIndex += 1
                    if indexLine.startswith('>'):
                        indexGuideCounter += 1
                        guideID = indexLine.split('>')[-1]
                        guideIDList.append(guideID)
            print('     There are {:,} guides in the index {}'.format(len(guideIDList), index_name))

            matchCounter = 0
            totalCounterSam = 0
            bamfileguideIDList = list() ### a list to hold all of the guides id found in sam file
            bam_file_path = working_directory + '/Bowtie2_Aligned_Filtered_files/' + file
            notFoundFile = bam_file_path.split('.bam')[0] + '_not_in_NGS_data.csv'
            reoccuranceFile = bam_file_path.split('.bam')[0] + '_occurrence.csv'
            with open(bam_file_path, 'r') as samReader:
                for bamLine in samReader:
                    totalCounterSam += 1
                    bamLine = bamLine.strip('\n')
                    bamLine = bamLine.split('\t')
                    if bamLine[2] == '*' or bamLine[0].startswith('@'):# *: non primary read and '@': header lines: should remove these
                        continue
                    else:
                        matchCounter += 1
                        bamfileguideIDList.append(bamLine[2])
            print('     Total number of primary reads in output sam file: {:,}'.format(matchCounter))

            not_in_list = list()  # a list to hold guide not found in the align data
            found_counter = 0  # count the number of guide found in the align data
            not_found_counter = 0
            temp_guideID_list = list()  # a temporary list to hold the keys from the counter dictionary
            for key,value in Counter(bamfileguideIDList).items():
                temp_guideID_list.append(key)
            with open(reoccuranceFile, 'w+') as re_occ_writer:
                for key,value in Counter(bamfileguideIDList).items():
                    line_to_write = key + ',' + str(value)
                    re_occ_writer.write(line_to_write + '\n')
            for guide in guideIDList:
                if guide in temp_guideID_list:
                    found_counter += 1
                else:
                    not_found_counter += 1
                    not_in_list.append(guide)

            with open(notFoundFile, 'w+') as not_found_writer:
                for element in not_in_list:
                    not_found_writer.write(element + '\n')
            coverage = found_counter / len(guideIDList) * 100
            uncoverage = not_found_counter/len(guideIDList) * 100

            print('     Number of guide found in NGS data: {:,}. Coverage rate: {:.2f}%'.format(found_counter, coverage))
            #print('     Number of guide not found in NGS data: {:,}. Uncoverage rate: {:.2f}%'.format(not_found_counter, uncoverage))


def step8_making_histogram(working_directory):
    '''Function to make the read counts vs frquency bar graph for the CRISPR analysis'''
    folder_location = working_directory + '/Bowtie2_Aligned_Filtered_files'
    for file in os.listdir(folder_location):
        if file.endswith('occurrence.csv') and 'unmatched' not in file:
            plot_title = file.split('_trimmed_loc')[0]
            file_path = folder_location + '/' + file
            count_list = list()  # a list to hold all of the count value
            bin_size = list()  # this is the temporary bin
            with open(file_path, 'r') as file_reader:
                for file_line in file_reader:
                    file_line = file_line.strip('\n')
                    file_line = file_line.split(',')
                    count_list.append(int(file_line[1]))

            a = np.array(count_list)
            x = 1
            step_bin = 66
            while x <= max(count_list):
                bin_size.append(x)
                x += step_bin
            bin_size.append(max(bin_size)+int(step_bin))  # temp_bin with another step_bin
            fig, ax = plt.subplots(figsize=(10, 7))
            ax.hist(a, bins=bin_size)
            plt.title(plot_title)
            plt.xlabel('Read Count', fontsize=15)
            plt.ylabel('Frequency', fontsize=15)
            histogramfilepath = folder_location + '/' + plot_title + '.png'
            plt.savefig(histogramfilepath)


def step9_combine_output_files(working_directory, index_name):
    '''Function to combine multiple occurrence output files after alignment'''
    foldernametocreate = 'Main_Outputs'
    MYDIR = working_directory + '/' + foldernametocreate
    CHECKFOLDER = os.path.isdir(MYDIR)
    if not CHECKFOLDER:
        os.makedirs(MYDIR)
        print('Created Folder', foldernametocreate)
    else:
        print('Folder', foldernametocreate, 'already exists!')

    folder_path = working_directory + '/Bowtie2_Aligned_Filtered_files/'
    output = working_directory + '/Main_Outputs/' + 'combined_output_' + index_name + '.csv'
    index_file_path = index_folder + index_name + '.fa'
    index_counter = 0  # count the number of guide in index file
    index = index_name.split('.')[0]
    index_list = list()  # a list to hold all guide from index
    with open(index_file_path, 'r') as indexreader:
        for indexline in indexreader:
            indexline = indexline.strip('\n')
            if indexline.startswith('>'):
                index_counter += 1
                indexline = indexline.replace('>','')
                index_list.append(indexline)
    print('Number of guide in index {}: {:,}'.format(index, index_counter))

    df = pd.DataFrame(index_list, columns=[index])

    for file in os.listdir(folder_path):
        if file.endswith('occurrence.csv') and 'unmatched' not in file and 'Bar' in file:
            file_dict = dict()  # a dictionary to hold key: guide name, value: occurrence
            file_path = folder_path + file
            reduced_file_name = file.split('_trimmed_')[0]
            with open(file_path) as reader:
                for line in reader:
                    line = line.strip('\n')
                    line = line.split(',')
                    file_dict[line[0]] = line[1]
            df[reduced_file_name] = df[index].map(file_dict)
    df = df.replace(np.nan, 0)
    df.to_csv(output, index=False)


def step10_prep_for_mageck(working_directory):
    '''Function to format count file to fit mageck software'''
    folder_path = working_directory + '/Main_Outputs/'
    for file in os.listdir(folder_path):
        if file.endswith('.csv') and file.startswith('combined') and 'mageck' not in file:
            inputfile = folder_path + file
            outputfile = inputfile.split('.csv')[0] + '_updated_for_mageck.csv'
            with open(outputfile, 'w+') as writer:
                with open(inputfile, 'r') as inputreader:
                    for line in inputreader:
                        line = line.strip('\n')
                        line = line.split(',')
                        if 'index' in line[0]:
                            line.insert(1, 'Gene')
                            linetowrite = ','.join(line)
                            writer.write(linetowrite + '\n')
                        else:
                            genename = line[0].split('|')[0]
                            line.insert(1, genename)
                            linetowrite = ','.join(line)
                            writer.write(linetowrite + '\n')


def step11_combine_output_deseq2(working_directory,index_name):
    '''Function to combine deseq2 R output and python output into 1 file'''
    foldernametocreate = 'Combined_with_DeSeq2'
    MYDIR = working_directory + '/Main_Outputs/' + foldernametocreate
    CHECKFOLDER = os.path.isdir(MYDIR)
    if not CHECKFOLDER:
        os.makedirs(MYDIR)
        print('Created Folder', foldernametocreate)
    else:
        print('Folder', foldernametocreate, 'already exists!')

    folder_path = working_directory + '/Main_Outputs'

    deseq2file = working_directory + '/Main_Outputs/deseq2_output.csv'
    python_output = working_directory + '/Main_Outputs/combined_output_'+ index_name + '.csv'
    combine_output = MYDIR + '/' + python_output.split('/')[-1].split('.csv')[0]+'_combined_with_deseq2.csv'

    python_file_line_counter = 0
    deseq2_file_line_counter = 0
    python_dict = dict()
    deseq2_dict = dict()
    with open(python_output, 'r') as python_reader:
        for python_line in python_reader:
            python_file_line_counter += 1
            python_line = python_line.strip('\n')
            python_line = python_line.split(',')
            python_line_value_list = python_line[1:]
            python_dict[python_line[0]] = ','.join(python_line_value_list)

    with open(deseq2file, 'r') as deseq2_reader:
        for deseq2_line in deseq2_reader:
            deseq2_file_line_counter += 1
            deseq2_line = deseq2_line.strip('\n')
            deseq2_line = deseq2_line.replace('"', '')
            deseq2_line = deseq2_line.split(',')
            if deseq2_line[0] == '':
                print(deseq2_line)
                deseq2_line[0] = index_name
                deseq2_line_value_list = deseq2_line[1:]
                deseq2_dict[deseq2_line[0]] = ','.join(deseq2_line_value_list)
            else:
                deseq2_line_value_list = deseq2_line[1:]
                deseq2_dict[deseq2_line[0]] = ','.join(deseq2_line_value_list)

    output_dict = defaultdict(list)

    for d in (python_dict, deseq2_dict):
        for key, value in d.items():
            output_dict[key].append(value)
    print(len(output_dict))

    with open(combine_output, 'w+') as output_writer:
        for k, v in output_dict.items():
            line_to_write = k + ',' + ','.join(v)
            output_writer.write(line_to_write + '\n')

    df = pd.read_csv(combine_output)
    df = df.replace(np.nan, 0)
    df['neg log10(pvalue)'] = np.negative(np.log10(df['pvalue']))
    df = df.replace(np.nan, 0)
    df = df.replace(np.inf, 0)
    df.to_csv(combine_output, index=False)


