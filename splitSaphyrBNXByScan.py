#!/usr/bin/env python
# Code Written by Jeffrey G. Reifenberger for BioNano Genomics
import re
import argparse
import sys

#function determines number of cohorts per scan in bnx file.
#if cohortStr is empty, ICS version should be 4 or greater. 8 cohorts per scan
#if cohortStr is not empty and contains the word 'Cohort', then ICS version is less than 4. 16 cohorts per scan
#if cohortStr is not empty and does not contain the word 'Cohort', then assume 8 cohorts per scan.
def determineCohortPerScan(cohortStr):
    if len(cohortStr) == 0: 
        return 8
    else:
        check = re.search('Cohort',cohortStr)
        if check:
            return 16
        else:
            return 8


#Function only reads bnx file header, i.e. lines that start with '#' in bnx file.
#Associates the runID in the header with the the correct scan number, last entry in lines that start with #Run Data
#Each molecule, line that starts with 0, contains a runID value as well.
#Function asigns a scan # for each runID in the header
#Later in code as each each molecule is read the runID of the molecule is associated with scan #.
#if ICS is version 4 or greater, there are 8 cohorts per scan. 
#if ICS is version less then 4, there are 16 cohorts per scan.
def readBNXHeader(filename):

    #read filename, i.e. bnx file.
    infile = open(filename,'rU')

    runIDHash = {}  #hash, key = runID value, value is list with scan number as value in list. [scanNum]
    numCohortPerScan = -1   #initialize at -1 so value can be determined in code based on header information.
    n = 0 #n tracks the number of Cohorts. When n = snumberCohortPerScan, reset to 0.
    scanNum = 1 #scanNum Tracks scan number of runID value.
    for line in infile: #iterate through each line in .bnx file.
        if line[0] == '#':  #only look at lines that are the header.
            if line[0:10] == '# Run Data':      #if line contains # Run Data,

                word = line.strip().split('\t') #split line based on tab, i.e. '\t'
                n+=1    #increment n

                if numCohortPerScan == -1:  #if numCohortPerScan = -1, determeCohortPerScan Value
                    numCohortPerScan = determineCohortPerScan(word[1])

                if n == numCohortPerScan:   #when n=numCohortPerScan, iterate scanNum, reset n
                    runIDHash[word[-1]] = [scanNum] #add runID and scanNum to runIDHash before iterate scanNum
                    scanNum+=1
                    n=0
                else:                       #else, add runID and scanNum to runIDHash
                    runIDHash[word[-1]] = [scanNum]
        else:       #only read header lines, #...., in bnx file. Once out of header return runIDHash.
            print '#### Code found '+str(scanNum-1)+' scans in bnx file.'
            print '#### Code found '+str(numCohortPerScan)+' number of cohorts per scan in bnx file.'
            return runIDHash    #Important: assume runId is proper order of data collection

#function prints how runID and scan numbers were assigned.
def printRunIDInformation(runIDHash,prefix):
    outfile = open(prefix+'_runID_to_scan.txt','w')

    outfile.write('RunID\tScanNumber\n')

    rkeys = runIDHash.keys()    #list of runID values.

    ikeys = []

    for r in rkeys: #rkeys is random list of runID values.
        ikeys.append(int(r))    #create list of runID values as integers so can sort.

    ikeys.sort()    #sort from runID=1 to runID=N

    for i in ikeys: #iterate through each r, print each r and scan number value.
        outfile.write(str(i)+'\t'+str(runIDHash[str(i)][0])+'\n')

    outfile.close()


#function determines if should split bnx file for each scan number or for give range of scan numbers
#returns scanHash that contains a string of potential Scan## names depending on how bnx file will be parsed.
def parseScanRange(sValue,runIDHash,strScanRange):
    #scanHash. key = scan number, value = Scan## (or if range of scan values selected, Scan##-Scan##)
    scanHash = {}

    #rkeys is list of runID values.
    rkeys = runIDHash.keys()

    if sValue == 1: #if sValue = 1, then generate bnx file for each scan number.
        for r in rkeys:     #iterate through each runID value, r
            s = str(runIDHash[r][0])    #runIDHash[r][0] is scan number

            if s in scanHash:   #if s in scanHash, donothing
                donothing = 1
            else:               #else, add Scan## to scanHash[s]
                scanHash[s] = 'Scan'+s.zfill(2)


    elif sValue == 0:       #if sValue = 0, then generate bnx file for range of scan numbers based on user input.
        rScan = strScanRange.split(',')     #for each ',' in strScanRange, split and place values into rScan.

        scanName = ''
        for r in rScan:
            scanName+=r+'_'     #scanName will contain str of value, e.g. if strScanRange = 1-10,15-20 then scanName = 1-10_15-20

        for r in rScan:         #for each range r (e.g. 1-10, 15-20, etc.) split value based on '-' to generate 1 or 2 numbers 
            sValue = r.split('-')

            if len(sValue) == 2:    #if split('-') results in 2 values in sValue, then interate between each value.
                for i in range(int(sValue[0]),int(sValue[1])+1):
                    if str(i) in scanHash:  #if i in scanHash, do nothing
                        donothing = 1
                    else:                   #else, add i to scanHash with value ScanRange_scanName (e.g. ScanRange_1-10)
                        scanHash[str(i)] = 'ScanRange_'+scanName
            elif len(sValue) == 1:  #if split('-') results in 1 values in sValue, then interate through the single value.
                for i in range(int(sValue[0]),int(sValue[0])+1):
                    if str(i) in scanHash:  #if i in scanHash, do nothing
                        donothing = 1
                    else:                   #else, add i to scanHash with value ScanRange_scanName (e.g. ScanRange_1)
                        scanHash[str(i)] = 'ScanRange_'+scanName

    return scanHash



def writeTitle(line,openFileHash,sValue):
    skeys = openFileHash.keys()     #list of scan numbers

    if sValue == 0:     #print single bnx file for a range of scan numbers. pick first skey value.
        openFileHash[skeys[0]].write(line)

    elif sValue == 1:     #print bnx file for each scan number.
        for s in skeys:
            openFileHash[s].write(line)


#generate hash of new bnx file names.
def generateFileNames(sValue,scanHash,runIDHash,prefix):
    #filenameHash: key=scan number, value is filename.
    filenameHash = {}

    rkeys = runIDHash.keys()

    #if sVavlue=0, then generate single bnx file for range of scan numbers
    #scanHash[s] = ScanRange_1-10_15-20 (for example)
    if sValue == 0:
        for r in rkeys:     #iterate through each runID
            s = str(runIDHash[r][0])    #runIDHash[r][0] is scan number

            if s in scanHash:
                if s in filenameHash:   #if s already in filenameHash, do nothing.
                    donothing = 1
                else:                   #else add s to filenameHash with new file name
                    filenameHash[s] = prefix+'_'+scanHash[s]+'filtered.bnx'

    #if sVavlue=1, then generate bnx file for each scan number
    #scanHash[s] = Scan01, Scan02, etc. for each s.
    elif sValue == 1:
        for r in rkeys:     #iterate through each runID
            s = str(runIDHash[r][0])    #runIDHash[r][0] is scan number

            if s in filenameHash:   #if s already in filenameHash, do nothing.
                donothing = 1
            else:                   #else add s to filenameHash with new file name
                filenameHash[s] = prefix+'_'+scanHash[s]+'.bnx'


    return filenameHash


def openFiles(sValue,filenameHash):
    #openFileHash: value is scan number, value is open file object.
    openFileHash = {}

    #if sValue = 0, then print range of scan numbers in a single bnx file.
    if sValue == 0:
        skeys = filenameHash.keys() #skeys is list of scan numbers to print in new bnx file.

        for i in range(0,len(skeys)):
            if i == 0:  #if i=0, open file, assign to fcode, place fcode in openFileHash[skeys[i]]
                fcode = open(filenameHash[skeys[0]],'w')    #open for write!
                openFileHash[skeys[i]] = fcode
            else:   
                openFileHash[skeys[i]] = fcode      #since only 1 file, assign fcode to each scan number, skeys[i]

    #if sValue = 1, then print a bnx file for each scan number.
    elif sValue == 1:
        skeys = filenameHash.keys()     #skeys is list of scan numbers to print in new bnx file.

        for s in skeys:
            openFileHash[s] = open(filenameHash[s],'w')     #open each bnx file per scan. for write!

    return openFileHash

#function reads bnx file and generates new bnx files based on user input.
def readBNXFile(sValue,scanHash,runIDHash,filename,prefix):
    infile = open(filename,'rU')    #open original bnx file. read only!

    #generate hash of new bnx file names based on user input.
    filenameHash = generateFileNames(sValue,scanHash,runIDHash,prefix)
    skeys = filenameHash.keys()     #fkeys is list of 

    openFileHash = openFiles(sValue,filenameHash)


    for line in infile: #iterate through each line in bnx file.
        if line[0] == '#':
            writeTitle(line,openFileHash,sValue)
        else:
            if line[0] == '0':              #check data from line 0
                keep = 0            #for each molecule, set keep=0, then check runID, scan number to determine if keep=1
                word = line.split() #split line 0 based on white space.

                runID = str(word[11])   #runID, which connects the molecule to the respective scan number
                snum = str(runIDHash[runID][0])    #scan number 

                if snum in openFileHash:    #if scan number in openFileHash, keep=1, write line in new bnx file.
                    keep = 1
                    openFileHash[snum].write(line)
                else:                       #else, scan number not in openFileHash, keep = 0, do not write line in new bnx file.
                    keep = 0

            elif line[0] == '1': #data from line 1, nick site locations, color 1
                if keep == 1:   #if keep = 1, then write line in new bnx file.
                    openFileHash[snum].write(line)

            elif line[0] == '2': #data from line 1, nick site locations, color 2
                if keep == 1:   #if keep = 1, then write line in new bnx file.
                    openFileHash[snum].write(line)

            elif line[0:4] == 'QX11': #data color 1 SNR values.
                if keep == 1:   #if keep = 1, then write line in new bnx file.
                    openFileHash[snum].write(line)

            elif line[0:4] == 'QX12': #data color 1 intensity values.
                if keep == 1:   #if keep = 1, then write line in new bnx file.
                    openFileHash[snum].write(line)

            elif line[0:4] == 'QX21': #data color 2 SNR values.
                if keep == 1:   #if keep = 1, then write line in new bnx file.
                    openFileHash[snum].write(line)

            elif line[0:4] == 'QX22': #data color 2 intensity values.
                if keep == 1:   #if keep = 1, then write line in new bnx file.
                    openFileHash[snum].write(line)


    skeys = openFileHash.keys()     #list of scan numbers in openFileHash.
    for s in skeys:
        openFileHash[s].close()     #close each file.


parser = argparse.ArgumentParser(description='Code filters bnx file based on scan number. Writes a new filtered BNX file based on input from user. Code should work for a 1 or 2 color bnx file.')
parser.add_argument("-b", "--bnxFile", help="full original bnx filename",type=str,default='')
parser.add_argument("-s", "--byScan", help="if -s 1, then print a bnx file for each scan #. if -s 0 print bnx file for range of scan numbers listed by user in -r. IMPORTANT: Must give values for -r if -s 0. default=0",type=int,default=0)
parser.add_argument("-r", "--scanRange", help="range of scan numbers to print new bnx file. For exmaple: 1-10,15-20 will print new bnx file with scans 1 through 10 and 15 through 20, but skip scans 11 through 14. default = '' ",type=str,default='')
parser.add_argument("-p", "--prefix", help="output bnx file prefix. default=some_great_data",type=str,default='some_great_data')
args = parser.parse_args()

#read bnx header, assign a san number to each runID.
runIDHash = readBNXHeader(args.bnxFile)

#print simple text file. Print each scan number, runID value.
printRunIDInformation(runIDHash,args.prefix)

#Check: If -s 0, code expects a value for -r. If not value give for -r exit code and print warning.
if args.byScan == 0:
    if len(args.scanRange) == 0:
        print
        print
        print '#### WARNING: You have selected to print bnx for range of scan numbers (-s 0), but have not given a range of scan numbers (-r '+args.scanRange+'). This is incompatible. If you select -s 0, then you must give a range of scan numbers (-r). For example: -r 1-10'
        sys.exit()

    elif len(args.scanRange) > 0:   #if values given to -r, continue
        #determine how to parse up the scan numbers based on the range of scan numbers give by user.
        scanHash = parseScanRange(0,runIDHash,args.scanRange)

#Check: If -s 1, then code will print a bnx file for each scan number. There should be no value for -r. 
elif args.byScan == 1:
    if len(args.scanRange) == 0:
        #parse bnx file for each scan number.
        scanHash = parseScanRange(1,runIDHash,args.scanRange)

    elif len(args.scanRange) > 0:   #if values listed for -r, then exit code and print warning.
        print
        print
        print '#### WARNING: You have selected to print a bnx file for each scan number (-s 1), but also given a range of scan numbers (-r '+args.scanRange+'). This is incompatible. If you select -s 1, then you can not select a range of scan numbers (-r) as well.'
        sys.exit()


print 'Reading, Filtering BNX File:\t',args.bnxFile
#read bnx file and generate new bnx file based on user input.
readBNXFile(args.byScan,scanHash,runIDHash,args.bnxFile,args.prefix)











