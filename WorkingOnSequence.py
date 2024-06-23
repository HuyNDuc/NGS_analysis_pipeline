# author: Huy Duc
# University of Colorado Denver - Anschutz Medical Campus
# Functional Genomics Facility

chromosome = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20',
          '21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37','38','MT','X']
databasePath = '/Users/huy.duc/OneDrive/Anschutz_Medical_Campus_Computational/Canfam3.1_database/'

import sqlite3
import re
import math
import gzip

def sequenceLookerReverse(loc,chrom):
    # chrom: name of the chromosome (1, 2, 3,...,X, MT)
    # location: an array of start number of the sequence
    sequenceFilePath = databasePath
    with open(sequenceFilePath + 'Full.sequences/sequence.' + chrom + '.txt', 'r') as read:
        for l in read:
            # counter = 1 #start at 1 for easy visualization
            # for loc in location:
            loc = int(loc)
            start = loc - 1 #20 mer + PAM
            end = loc + 22 #20 mer + PAM
            seq = l[start:end]
            seq30mer = l[start - 3:end + 4] #30 mer: 4bps + 20 mer + PAM + 3bps
            # print(counter,'::',l[start:end])
            # counter += 1
            result = [get_reverse_compliment(seq), get_reverse_compliment(seq30mer)]  # return[0]: 20mer + PAM, [1]: 30mer
    return result

# This function is for looking for 30_mers
def sequenceLooker2(loc,chrom):
    with open(databasePath+'Full.sequences/sequence.' + chrom + '.txt', 'r') as read:
        for l in read:
            # counter = 1 #start at 1 for easy visualization
            # for loc in location:
            start = loc
            end = loc + 30
            seq = l[start:end]

    return seq

def get_reverse_compliment(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    letters = list(s[::-1])
    letters = [basecomp.get(base, 'N') for base in letters]
    return ''.join(letters)

def get_revese(s):
    letters = list(s[::-1])
    return ''.join(letters)

def get_compliment(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A'}
    letters = list(s[::1])
    letters = [basecomp.get(base, 'N') for base in letters]
    return ''.join(letters)

def get_RNA(s):
    basecomp = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A'}
    letters = list(s[::1])
    letters = [basecomp.get(base, 'N') for base in letters]
    return ''.join(letters)


def dnaSearcherExactMatch(): #function to find exact match, only work for forward sequece

    for i in chromosome:
        print('reading chromosome ' + i +' ...')
        dbPath = "SequenceDataSorted\massiveSequenceForChrom." + i + ".db"
        connection = sqlite3.connect(dbPath)
        cursorF = connection.execute("SELECT g20_mer, start_20g, end_20g  from F_Pam")
        cursorR = connection.execute("SELECT g20_mer, start_20g, end_20g  from R_Pam")
        seqStr = "AAGAATCGAATCGAATCGAA"

        for data in cursorF:
            if data[0] == seqStr:
                print('     found in chromosome ' + i +', '+ str(data[1]) +':' +str(data[2]) + ' forward sequence')
        for data in cursorR:
            if get_revese(data[0]) == seqStr:
                print('     found in chromosome ' + i +', '+ str(data[1]) +':' +str(data[2]) + ' reverse sequence')


def getDataFromCDS(gene):
    # function to extract information of given gene from database generated from CDSf
    dbPath = 'CDSdatabase.db'
    count = 0
    connection = sqlite3.connect(dbPath)
    cursor = connection.execute("SELECT gene_name, exon_number, sequence, rev_sequence  from massiveExonDataV4")
    for data in cursor:
        if data[0].upper() == gene.upper():
            count = count + 1
            print(data[0],':','exon ', data[1],':')
            print('Exon length: ' + str(len(data[2]))+' bps')
            print('From forward sequence:' + data[2])
            get20Mers(data[2])
            print('From reversed sequence:' + data[3])
            revSeq = get_revese(data[3])
            get20Mers(revSeq)
        else:
            continue
    connection.commit()
    connection.close()
    if count == 0:
        print('Gene ' + gene + ' was not found in the database!')
    else:
        print('Gene ' + gene + ' has ' + str(count) + ' exons')


#this function is not working well, using bowtie instead
def aligning(a,b): #function to align two sequence with same length, from biopython
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(a, b)
    score = aligner.score(a, b)
    result = [score, alignments[0]] #return an array of two results, score and align
    return result

#function to test aligning pattern and genome
def testAligning():
    testChrom = ['1']
    seq1 = 'TTGAGAGCTGAGGACAGGAC'
    matchChrom4 = []
    matchChrom5 = []
    matchChrom6 = []
    match = []
    for t in testChrom:
        testCounterFoward = 0
        testCounterReverse = 0
        FLoc = []
        RLoc = []
        dbPath = "SequenceDataSorted\massiveSequenceForChrom."+t+".db"
        connection = sqlite3.connect(dbPath)
        foward = connection.execute('SELECT g20_mer, start_20g, end_20g FROM F_PAM')
        reverse = connection.execute('SELECT g20_mer, start_20g, end_20g FROM R_PAM')
        print('reading forward strand for chromosome', t)
        for f in foward:
            aResult = aligning(seq1, f[0])
            if aResult[0] >= 18: #20 would be a perfect match
                print(aResult[0])
                #print(aResult[1])
                testCounterFoward += 1
                print('location forward:', f[1], ':', f[2]) #location of the match
                strLocF = str(str(f[1])+':'+str(f[2]))
                FLoc.append(strLocF)
                Flocstr = ', '.join(map(str,FLoc))
        print('reading reversed strand for chromosome', t)
        for r in reverse:
            rev_seq = get_revese(r[0]) # because it is reversed strand
            rResult = aligning(seq1, rev_seq)
            if rResult[0] >= 18:
                print(rResult[0])
                #print(rResult[1])
                testCounterReverse += 1
                print('location reverse:', f[1], ':', f[2]) #location of the match
                strLocR = str(str(f[1])+':'+str(f[2]))
                RLoc.append(strLocR)
                Rlocstr = ', '.join(map(str,RLoc))
        if len(FLoc) == 0:
            print('No match found in chromosome',t,'foward strand')
        else:
            print(testCounterFoward, 'match(es) found in chromosome',t,'at',Flocstr)
        if len(RLoc) == 0:
            print('No match found in chromosome',t,'reverse strand')
        else:
            print(testCounterReverse, 'match(es) found in chromosome', t, 'reverese strand at',Rlocstr)

        connection.commit()
        connection.close()

def exactMatch(pat,chrom): #not as good as bowtie, use bowtie, this method obmits some of the correct results
    with open(databasePath+'Full.sequences/sequence.'+chrom+'.txt', 'r') as read:
        for e in read:
            resultF = [m.start() for m in re.finditer(pat, e)]
            eR = get_reverse_compliment(e)
            resultR = [m.start() for m in re.finditer(pat, eR)]
            print(resultF, len(resultF))
            print(resultR, len(resultR))
            print('found', len(resultF), 'matches in forward sequence of chromosome',chrom)
            print('found', len(resultR), 'matches in reverse sequence of chromosome',chrom)

def generatefqFileFromCDSForward():
    #create fq file for guides in forward strand, with gene name and exon number
    dbPath = 'CDSdatabase.db'
    count = 0
    connection = sqlite3.connect(dbPath)
    cursor = connection.execute("SELECT sequence, gene_name, exon_number, ensemblID  from massiveExonDataV4")
    fqFile = open('CanFam3.1.96GuideForward.fq', 'w+')
    for data in cursor:
        #a simplified version of get20mers
        pam = 'CGG|AGG|TGG|GGG'
        findPam = [m.start() for m in re.finditer(pam, data[0])]
        endLim = len(data[0])
        for i in findPam:
            start = i - 20
            end = i
            if start >= 0 and end <= endLim:
                s_20mers = data[0][start:end]
                testPAM = data[0][start:end+3]
                fqFile.write('@CF3.1.96-'+str(count)+'::'+data[3]+'::'+str(data[1])+'::'+str(data[2])+'::AGG'+'\n')
                fqFile.write(s_20mers+'AGG'+'\n')
                fqFile.write('+'+'\n')
                fqFile.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                count += 1
                fqFile.write('@CF3.1.96-'+str(count)+'::'+data[3]+'::'+str(data[1])+'::'+str(data[2])+'::GGG'+'\n')
                fqFile.write(s_20mers + 'GGG'+'\n')
                fqFile.write('+'+'\n')
                fqFile.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                count += 1
                fqFile.write('@CF3.1.96-'+str(count)+'::'+data[3]+'::'+str(data[1])+'::'+str(data[2])+'::TGG'+'\n')
                fqFile.write(s_20mers + 'TGG'+'\n')
                fqFile.write('+'+'\n')
                fqFile.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                count += 1
                fqFile.write('@CF3.1.96-'+str(count)+'::'+data[3]+'::'+str(data[1])+'::'+str(data[2])+'::CGG'+'\n')
                fqFile.write(s_20mers + 'CGG'+'\n')
                fqFile.write('+'+'\n')
                fqFile.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                count += 1
                print(s_20mers, len(s_20mers), testPAM, count)

    fqFile.close()
    connection.commit()
    connection.close()

def generatefqFileFromCDSReverse():
    dbPath = 'CDSdatabase.db'
    count = 0
    connection = sqlite3.connect(dbPath)
    cursor = connection.execute("SELECT sequence  from massiveExonDataV4")
    fqFileR = open('CanFamGuideReverse.fq', 'w+')
    for data in cursor:
        #a simplified version of get20mers
        revPam = 'GGC|GGA|GGT|GGG'
        compSeq = get_compliment(data[0])
        findRevPam = [t.start() for t in re.finditer(revPam, compSeq)]
        endLim = len(compSeq)
        for i in findRevPam:
            start = i + 3
            end = i + 23
            if start >= 0 and end <= endLim:
                s_20mers = compSeq[start:end]
                testPAM = compSeq[start-3:end]
                fqFileR.write('@CF3.1.96-'+str(count)+':::GGA, reverse'+'\n')
                fqFileR.write('GGA' + s_20mers+'\n')
                fqFileR.write('-'+'\n')
                fqFileR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                fqFileR.write('@CF3.1.96-'+str(count)+':::GGG, reverse'+'\n')
                fqFileR.write('GGG' + s_20mers+'\n')
                fqFileR.write('-'+'\n')
                fqFileR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                fqFileR.write('@CF3.1.96-'+str(count)+':::GGT, reverse'+'\n')
                fqFileR.write('GGT' + s_20mers+'\n')
                fqFileR.write('-'+'\n')
                fqFileR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                fqFileR.write('@CF3.1.96-'+str(count)+':::GGC, reverse'+'\n')
                fqFileR.write('GGC' + s_20mers+'\n')
                fqFileR.write('-'+'\n')
                fqFileR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                print(s_20mers, len(s_20mers), testPAM, count)
                count += 1
    fqFileR.close()
    connection.commit()
    connection.close()


def generatefqFileFromCDSForwardAndReverse():
    dbPath = 'CDSdatabase.db'
    count = 0
    connection = sqlite3.connect(dbPath)
    cursor = connection.execute("SELECT sequence  from massiveExonDataV4")
    fqFileFR = open('CanFamGuideForwardAndReverse.fq', 'w+')
    for data in cursor:
        #a simplified version of get20mers
        pam = 'CGG|AGG|TGG|GGG'
        findPam = [m.start() for m in re.finditer(pam, data[0])]
        revPam = 'GGC|GGA|GGT|GGG'
        compSeq = get_compliment(data[0])
        findRevPam = [t.start() for t in re.finditer(revPam, compSeq)]
        endLim = len(compSeq)
        for i in findPam:
            start = i - 20
            end = i
            if start >= 0 and end <= endLim:
                s_20mers = data[0][start:end]
                testPAM = data[0][start:end+3]
                fqFileFR.write('@CF3.1.96-'+str(count)+':::AGG'+'\n')
                fqFileFR.write(s_20mers+'AGG'+'\n')
                fqFileFR.write('+'+'\n')
                fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                fqFileFR.write('@CF3.1.96-'+str(count)+':::GGG'+'\n')
                fqFileFR.write(s_20mers + 'GGG'+'\n')
                fqFileFR.write('+'+'\n')
                fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                fqFileFR.write('@CF3.1.96-'+str(count)+':::TGG'+'\n')
                fqFileFR.write(s_20mers + 'TGG'+'\n')
                fqFileFR.write('+'+'\n')
                fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                fqFileFR.write('@CF3.1.96-'+str(count)+':::CGG'+'\n')
                fqFileFR.write(s_20mers + 'CGG'+'\n')
                fqFileFR.write('+'+'\n')
                fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                print(s_20mers, len(s_20mers), testPAM, count)
                count += 1
        for i in findRevPam:
            start = i + 3
            end = i + 23
            if start >= 0 and end <= endLim:
                s_20mersR = compSeq[start:end]
                testPAMR = compSeq[start-3:end]
                fqFileFR.write('@CF3.1.96-'+str(count)+':::GGA, reverse'+'\n')
                fqFileFR.write('GGA' + s_20mersR+'\n')
                fqFileFR.write('-'+'\n')
                fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                fqFileFR.write('@CF3.1.96-'+str(count)+':::GGG, reverse'+'\n')
                fqFileFR.write('GGG' + s_20mersR+'\n')
                fqFileFR.write('-'+'\n')
                fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                fqFileFR.write('@CF3.1.96-'+str(count)+':::GGT, reverse'+'\n')
                fqFileFR.write('GGT' + s_20mersR+'\n')
                fqFileFR.write('-'+'\n')
                fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                fqFileFR.write('@CF3.1.96-'+str(count)+':::GGC, reverse'+'\n')
                fqFileFR.write('GGC' + s_20mersR+'\n')
                fqFileFR.write('-'+'\n')
                fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII'+'\n')
                print(s_20mersR, len(s_20mersR), testPAMR, count, ' ::: reverse sequence')
                count += 1
    fqFileFR.close()
    connection.commit()
    connection.close()


def readMapFie(file): #function to read the map file generated by bowtie
    with open(file, 'r') as f:
        count = 1
        for line in f:
            lineArray = line.split('\t')
            count += 1
            #result array:[0]: ID, [1]:strand, [2]:ref, [3]:location, [4]:sequence, [5]:quality, [6]:# of matches, [7]:new line
            print(count)

# def mapFileAnalyzer():
#     # This function is not working with mac, need to fix path information
#     with open('pak3_3mis.map', 'r') as iMap:
#         idDictionary ={}
#         for line in iMap:
#             line = line.split('\t')
#             #0:ID, [1]:strand, [2]:chrom, [3]:start locatioin, [4]:sequence, [5]:quality, [6]:??, [7] mismatch
#             idDictionary.update({line[0]:[line[3]]})
#     with open('pak3_3mis.map','r') as iMap:
#         for line in iMap:
#             line = line.split(('\t'))
#             for k,v in idDictionary.items():
#                 if line[0] == k:
#                     if line[3] in v: #so that the dictionary does not add in the already exist location
#                         continue
#                     else:
#                         idDictionary[k].append(line[3])
#                 else:
#                     continue
#     print('finish setting up dictionary for pak3 gene')
#
#     dbPath = dataPath + 'CDSdatabaseCanFam3.1.96.db'
#     connection = sqlite3.connect(dbPath)
#     cursor = connection.execute("SELECT chromosome, gene_name, CDS_start, CDS_end from massiveExonDataV4")
#     locArray = []
#     for cur in cursor:
#         if cur[0] == 'X':
#             interArray = [cur[1],int(cur[2]),int(cur[3])]
#             #[0]: gene name, [1]: start, [2]: end
#             locArray.append(interArray)
#     connection.commit()
#     connection.close()
#     print('finish getting location data for CDS for chromosome X')
#
#     idDictionary = idDictionary.items()
#     for item in idDictionary:
#         print(item[0],item[1])
#         for i in item[1]:
#             i = int(i)
#             for loc in locArray:
#                 if i >= loc[1] and i <= loc[2]:
#                     print(i,'::this guide hit gene',loc[0],'between',loc[1],'and',loc[2])

def extractSampleData(): #extract some location data for testing
    #test: CF3.1.96-66:::CGG on PAK3, chromosome X,
    chr = 'X'
    #locDict = {}
    with open('pak3_3mis.map', 'r') as iMap:
        for line in iMap:
            line = line.split('\t')
            # 0:ID, [1]:strand, [2]:chrom, [3]:start location, [4]:sequence, [5]:quality, [6]:??, [7] mismatch
            if line[0] == 'CF3.1.96-66:::CGG':
                chrSeq = sequenceLooker(line[3],chr)
                #locDict.update({line[3]:chrSeq})

                print(line[3],':',line[1],':',line[4],':',line[2],':',chrSeq,line[7])

def getData(gene):
    dbPath = 'CanFam3.1exondata.v2.db'
    count = 0
    connection = sqlite3.connect(dbPath)
    cursor = connection.execute("SELECT gene_name, strand, sequence  from massiveExonData")
    for data in cursor:
        if data[0].upper() == gene.upper():
            if data[1] == '+':
                count = count + 1
                print(data[0] + ' : ' + data[1] + ' : ' + data[2])
                print('Exon length: ' + str(len(data[2]))+ ' bps')
                get20Mers(data[2])
            else:
                continue
            # elif data[1] == '-': #need to check this to make sure
            #     count = count + 1
            #     newseq = get_reverse_compliment(data[2])
            #     print(data[0] + ' : ' + data[1] + ' : ' + 'reversed: ' + newseq)
            #     getgRNA(newseq)
        else:
            continue
    connection.close()
    if count == 0:
        print('Gene ' + gene + ' was not found in the database!')
    else:
        print('gene ' + gene + ' has ' + str(count) + ' exons')


def getDataForAll():
    dbPath = 'CanFam3.1exondata.v2.db'
    connection = sqlite3.connect(dbPath)
    cursor = connection.execute("SELECT gene_name, strand, sequence  from massiveExonData")
    for data in cursor:
        print(data[0] + ' : ' + data[2])
        get20Mers(data[1])
    connection.close()


def findMatches(seq):
    for i in chromosome:
        x = ''
        with open(databasePath +'Full.sequences/sequence.'+i+'.txt', 'r') as read:
            for line in read:
                if seq in line:
                    x = i
                else:
                    x = ''
    return x

#this function calculate score using Rule Set 1, Doench 2014
def calc_score(s):
    s_list = list(s)
    s_20mer = s[4:24]
    nuc_hash = {'A':0, 'T':1, 'C':2, 'G':3}
    score = 0.597636154
    gc = s_20mer.count('G')+s_20mer.count('C')
    gc_low = -0.202625894
    gc_high = -0.166587752
    if gc < 10:
        gc_val = abs(gc-10)
        score = score+(gc_val*gc_low)
    elif gc > 10:
        gc_val = gc-10
        score = score+(gc_val*gc_high)
    #rows[1-30]cols['ATCG']
    sing_nuc_hash = {'G2':-0.275377128,'A3':-0.323887456,'C3':0.172128871,'C4':-0.100666209,'C5':-0.20180294, \
                    'G5':0.245956633,'A6':0.036440041,'C6':0.098376835,'C7':-0.741181291,\
                    'G7':-0.393264397,'A12':-0.466099015,'A15':0.085376945,'C15':-0.013813972,\
                    'A16':0.272620512,'C16':-0.119022648,'T16':-0.285944222,'A17':0.097454592,\
                    'G17':-0.17554617,'C18':-0.345795451,'G18':-0.678096426,'A19':0.22508903,\
                    'C19':-0.507794051,'G20':-0.417373597,'T20':-0.054306959,'G21':0.379899366,\
                    'T21':-0.090712644,'C22':0.057823319,'T22':-0.530567296,'T23':-0.877007428,\
                    'C24':-0.876235846,'G24':0.278916259,'T24':-0.403102218,'A25':-0.077300704,\
                    'C25':0.287935617,'T25':-0.221637217,'G28':-0.689016682,'T28':0.117877577,\
                    'C29':-0.160445304,'G30':0.386342585}
    #score_mat = np.matrix('0 0 0 0;0 0 0 -0.275377128;-0.323887456 0 0.172128871 0;0 0 -0.100666209 0;0 0 -0.20180294 0.245956633;0.036440041 0 0.098376835 0;0 0 -0.741181291 -0.393264397;0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0;-0.466099015 0 0 0;0 0 0 0;0 0 0 0;0.085376945 0 -0.013813972 0;0.272620512 -0.285944222 -0.119022648 0;0.097454592 0 0 -0.17554617;0 0 -0.345795451 -0.678096426;0.22508903 0 -0.507794051 0;0 -0.054306959 0 -0.417373597;0 -0.090712644 0 0.379899366;0 -0.530567296 0.057823319 0;0 -0.877007428 0 0;0 -0.403102218 -0.876235846 0.278916259;-0.077300704 -0.221637217 0.287935617 0;0 0 0 0;0 0 0 0;0 0.117877577 0 -0.689016682;0 0 -0.160445304 0;0 0 0 0.386342585')
    dinuc_hash = {'GT2':-0.625778696,'GC5':0.300043317,'AA6':-0.834836245,'TA6':0.760627772,'GG7':-0.490816749,
                  'GG12':-1.516907439,'TA12':0.7092612,'TC12':0.496298609,'TT12':-0.586873894,'GG13':-0.334563735,
                  'GA14':0.76384993,'GC14':-0.53702517,'TG17':-0.798146133,'GG19':-0.66680873,'TC19':0.353183252,
                  'CC20':0.748072092,'TG20':-0.367266772,'AC21':0.568209132,'CG21':0.329072074,'GA21':-0.836456755,
                  'GG21':-0.782207584,'TC22':-1.029692957,'CG23':0.856197823,'CT23':-0.463207679,'AA24':-0.579492389,
                  'AG24':0.649075537,'AG25':-0.077300704,'CG25':0.287935617,'TG25':-0.221637217,'GT27':0.117877577,
                  'GG29':-0.697740024}
    for i, nuc in enumerate(s_list):
        key = nuc+str(i+1)
        if key in sing_nuc_hash:
            nuc_score = sing_nuc_hash[key]
        else:
            nuc_score = 0
        #nuc_score = score_mat[i,nuc_hash[nuc]]
        #print(key + ' : ' + str(nuc_score))
        score = score+nuc_score
        if i < 29:
            dinuc = nuc+s[i+1]+str(i+1)
            if dinuc in dinuc_hash.keys():
                score = score+dinuc_hash[dinuc]
    partial_score = math.e**-score
    final_score = 1/(1+partial_score)
    return final_score


def countReoccurance(seq):
    counttotal = 0
    array = []
    resultdic = {}
    for i in chromosome:
        with open(databasePath+'Full.sequences/sequence.' + i + '.txt', 'r') as read:
            for line in read:
                count = line.count(seq)
                # print(i + ' : ' + str(count))
                if count > 0:
                    array.append(i)
            counttotal = counttotal + count
    resultdic.update({'totalcount' : counttotal})
    resultdic.update({'chrom' : array})
    return resultdic


def get20Mers(seq):
    pam = 'CGG|AGG|TGG|GGG'
    #guidefound = 0
    findPam = [m.start() for m in re.finditer(pam, seq)]
    endLim = len(seq)
    for i in findPam:
        start = i-24
        end = i + 6
        if start >= 0 and end <= endLim:
            s_30mers = seq[start:end]
            s_20mers = s_30mers[4:24]
            pam = s_30mers[24:27]
            score = calc_score(s_30mers)

            print(s_30mers + ", 20_mer sequence: " + s_20mers + ', PAM: ' + pam + ', Doench Score: ' + str(score))
            # if score >= 0.6:
            #     score = "%.2f" % score
            #     guidefound = guidefound + 1
            #     #print(s_30mers + ", guide sequence: " + gRNAseq + ', PAM: ' + pam + ', On Target Efficiency: ' + score)
            #     appearance = countReoccurance(gRNAseq)
            #     if appearance['totalcount'] == 1:
            #         counter = counter + 1
            #         print('\t' + str(counter) + ': ' + s_30mers + ' [' + str(start) +':'+str(end) +'] :' + gRNAseq + ' : PAM:' + pam + ' score = ' + str(score))
            #     else:
            #         print('\t' + gRNAseq + ' : matched somewhere in chromosome ' + str(appearance['chrom']))
            # else:
            #     continue
        else:
            continue
    # if guidefound == 0:
    #     print('\t' + 'No suitable guide was found for this exon')
    # elif guidefound == 1:
    #     print('\t' + str(guidefound) + ' suitable guide was found for this exon')
    # else:
    #     print('\t' + str(guidefound) + ' suitable guides were found for this exon')

def sequenceLooker(loc, chrom): #look for sequence with location within chromosome
    #chrom: name of the chromosome (1, 2, 3,...,X, MT)
    #location: an array of start number of the sequence
    sequenceFilePath = databasePath
    with open(databasePath+'/Full.sequences/sequence.' + chrom + '.txt', 'r') as read:
        for l in read:
            # counter = 1 #start at 1 for easy visualization
            # for loc in location:
            loc = int(loc)
            start = loc - 1
            end = loc + 22
            seq = l[start:end]
            seq30mer = l[start - 4:end + 3]
            # print(counter,'::',l[start:end])
            # counter += 1
            result = [seq, seq30mer] #return[0]: 20mer + PAM, [1]: 30mer
    return result


def getGeneFASTQ():
    geneToGet = 'PAK3'
    dbPath = 'CDSdatabase.db'
    count = 0
    connection = sqlite3.connect(dbPath)
    cursor = connection.execute("SELECT sequence, gene_name  from massiveExonDataV4")
    outputFASTQ = 'testFilePAK3.fq'
    fqFileFR = open(outputFASTQ, 'w+')
    for data in cursor:
        if data[1] == geneToGet:
            # a simplified version of get20mers
            pam = 'CGG|AGG|TGG|GGG'
            findPam = [m.start() for m in re.finditer(pam, data[0])]
            revPam = 'GGC|GGA|GGT|GGG'
            compSeq = get_compliment(data[0])
            findRevPam = [t.start() for t in re.finditer(revPam, compSeq)]
            endLim = len(compSeq)
            for i in findPam:
                start = i - 20
                end = i
                if start >= 0 and end <= endLim:
                    s_20mers = data[0][start:end]
                    testPAM = data[0][start:end + 3]
                    fqFileFR.write('@CF3.1.96-' + str(count) + ':::AGG' + '\n')
                    fqFileFR.write(s_20mers + 'AGG' + '\n')
                    fqFileFR.write('+' + '\n')
                    fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII' + '\n')
                    fqFileFR.write('@CF3.1.96-' + str(count) + ':::GGG' + '\n')
                    fqFileFR.write(s_20mers + 'GGG' + '\n')
                    fqFileFR.write('+' + '\n')
                    fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII' + '\n')
                    fqFileFR.write('@CF3.1.96-' + str(count) + ':::TGG' + '\n')
                    fqFileFR.write(s_20mers + 'TGG' + '\n')
                    fqFileFR.write('+' + '\n')
                    fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII' + '\n')
                    fqFileFR.write('@CF3.1.96-' + str(count) + ':::CGG' + '\n')
                    fqFileFR.write(s_20mers + 'CGG' + '\n')
                    fqFileFR.write('+' + '\n')
                    fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII' + '\n')
                    print(s_20mers, len(s_20mers), testPAM, count)
                    count += 1

#function to search for exact match
def queryExactMatch():
    fileToTest = 'MT#003A400-6000_within.csv'
    sequenceResult = 'MT400-6000seq.fa'
    #sequenceFile = "C:\Users\huydu\OneDrive\Anschutz Medical Campus\CRISPRlaptop\CanFam3.1.Chomosomes\Canis_familiaris.CanFam3.1.dna.chromosome.MT.txt"
    counter = 0
    with open(fileToTest,'r') as f:
        with open(sequenceResult, 'w+') as writer:
            for l in f:
                l = l.split(',')
                if l[0] == "MT":
                        print(sequenceLooker(l[1],'MT')[1])
                        # writer.write('@r'+str(counter)+'\n')
                        # writer.write(sequenceLooker(l[1],'MT')[1]+'\n')
                        # writer.write(l[6]+'\n')
                        # writer.write('IIIIIIIIIIIIIIIIIIIIIIIIIIIIII'+'\n')

                        writer.write('@r'+str(counter)+'\n')
                        writer.write(l[3]+"\n")
                        writer.write('+'+'\n')
                        writer.write("IIIIIIIIIIIIIIIIIIII"+'\n')
                        counter += 1

#function to generate fa file with all chromosome that could be used for guidescan
def generateFAfilesForAllChrom():
    with open('Canis_familiaris.CanFam3.1.dna.chromsome.all.fa','w+') as writer:
        for chr in chromosome:
            generalFilePath = 'C:\\Users\\huydu\\OneDrive\\Anschutz Medical Campus\\CRISPRlaptop\\CanFam3.1.Chomosomes\\gz files\\'
            filePath = generalFilePath + 'Canis_familiaris.CanFam3.1.dna.chromosome.'+chr+'.fa.gz'
            with gzip.open(filePath,'r') as fileContent:
                for line in fileContent:
                    line = str(line)
                    line = line.strip('b')
                    line = line.strip('\'')
                    line = line.strip('\\n')
                    writer.write(line+'\n')
            print(chr)

def checkChromLength(): #function to check the length of all chromosomes, write to a txt file
    with open(databasePath+'CanFam3.1.96-chromosomeslength.txt', 'w') as write:
        for i in chromosome:
            with open(databasePath+'Full.sequences/sequence.' + i + '.txt', 'r') as read:
                for l in read:
                    length = len(l)
                    length = f'{length:,}'
                    write.write(i +':::'+length+'\n')
                    print(i,'::',length)

#this function creates fastq file for each chromosome. used for testing purposes
def getFASTQFileForChrom(chroms):
    # arg chroms must be an array:
    for chrom in chroms:
        dbPath = 'CDSdatabase.db'
        count = 0
        connection = sqlite3.connect(dbPath)
        cursor = connection.execute("SELECT sequence, chromosome  from massiveExonDataV4")
        outputFASTQ = 'CanFamGuideForChrome.'+chrom+'.fq'
        fqFileFR = open(outputFASTQ, 'w+')
        for data in cursor:
            if data[1] == chrom:
                # a simplified version of get20mers
                pam = 'CGG|AGG|TGG|GGG'
                findPam = [m.start() for m in re.finditer(pam, data[0])]
                revPam = 'GGC|GGA|GGT|GGG'
                compSeq = get_compliment(data[0])
                findRevPam = [t.start() for t in re.finditer(revPam, compSeq)]
                endLim = len(compSeq)
                for i in findPam:
                    start = i - 20
                    end = i
                    if start >= 0 and end <= endLim:
                        s_20mers = data[0][start:end]
                        testPAM = data[0][start:end + 3]
                        fqFileFR.write('@CF3.1.96-' + str(count) + ':::AGG' + '\n')
                        fqFileFR.write(s_20mers + 'AGG' + '\n')
                        fqFileFR.write('+' + '\n')
                        fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII' + '\n')
                        fqFileFR.write('@CF3.1.96-' + str(count) + ':::GGG' + '\n')
                        fqFileFR.write(s_20mers + 'GGG' + '\n')
                        fqFileFR.write('+' + '\n')
                        fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII' + '\n')
                        fqFileFR.write('@CF3.1.96-' + str(count) + ':::TGG' + '\n')
                        fqFileFR.write(s_20mers + 'TGG' + '\n')
                        fqFileFR.write('+' + '\n')
                        fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII' + '\n')
                        fqFileFR.write('@CF3.1.96-' + str(count) + ':::CGG' + '\n')
                        fqFileFR.write(s_20mers + 'CGG' + '\n')
                        fqFileFR.write('+' + '\n')
                        fqFileFR.write('IIIIIIIIIIIIIIIIIIIIIII' + '\n')
                        print(s_20mers, len(s_20mers), testPAM, count)
                        count += 1
