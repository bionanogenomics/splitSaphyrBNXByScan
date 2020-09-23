#!/usr/bin/env python
# Code Written by Jeffrey G. Reifenberger for BioNano Genomics
# Updated 1.1 version, date 2020_03_11 that will better handle different bnx header formats.
# Updated 1.2 version, date 2020_09_23 will add the correct number of molecules for new bnx file in the header.
import re
import argparse
import sys
import numpy as np

#function will iterate through every molecule in .bnx file.
#It will search for the maximum colID for each value of runID.
#the maximum colID value for each runID informs the code on how to split the
#bnx file by scan.
def findNumberColumnsPerBank(filename):
    #open .bnx file, filename
    infile = open(filename,'rU')

    #list of runID values 
    runIDList = []

    #molNumPerRunID: keys is runID, value is number of molecules in runID
    molNumPerRunID = {}

    #runIDMaxColHash: key runID, value maximum colID for given runID.
    runIDMaxColHash = {}

    #maximum colID for entire .bnx file.
    maxColID = 0

    #iterate through each line in the .bnx file.
    for line in infile:
        #if the line begins with '# Run Data', pull out runID values.
        if line[0:10] == '# Run Data':
            #split line based on tab. add to word list.
            word = (line.strip()).split('\t')

            #add runID value, word[-1], to runIDList
            runIDList.append(int(word[-1]))
            #add runID, word[-1], to runIDMaxColHash, set value to 0. 
            runIDMaxColHash[str(word[-1])] = 0

        #read each line that begins with 0, molecule information.
        if line[0] == '0':
            #split line by whitespace. add values to word list.
            word = line.split()

            #pull runID and colID from 0 line.
            runID = str(word[11])
            colID = int(word[12])

            if runID in molNumPerRunID:
                molNumPerRunID[runID] += 1
            else:
                molNumPerRunID[runID] = 1


            #Check if colID is largest value observed for runID via runIDMaxColHash[runID].
            #if largest value, replace runIDMaxColHash[runID] with current colID value.
            if colID > runIDMaxColHash[runID]:
                runIDMaxColHash[runID] = colID

            #check if colID is largest value observed for all molecules. Store value in maxColID
            if colID > maxColID:
                maxColID = colID

    infile.close()

    #return maxColID,runIDList,runIDMaxColHash,molNumPerRunID
    return maxColID,runIDList,runIDMaxColHash,molNumPerRunID

#based on maximum colID per runID compute the number of runIDs per scan.
def determineRunIDPerScan(maxColID,runIDList,runIDMaxColHash):
    #list of computed number runIDs per bank on how many times the repetition of the maxColID value.
    runIDPerBank = []

    #number of banks in a flow cell. Hard coded to be 4.0 for now.
    numBanks = 4.0

    #a split FOV saphyr will have 137 imaging columns per bank --- sometimes the 137 columns are split to be 34, 34, 34, and 35 columnIDs per bank.
    #a full FOV saphyr will have 69 imaging columns per bank.
    #allow for some wiggle room regarding max number of columns per bank detected if chip is bad or loading is low.
    if maxColID > 69:   #this case is for split FOV saphyr with 137 image columns per bank.
        maxColID = 137
    elif maxColID <= 69 and maxColID > 35:      #this case is for full FOV saphyr with 69 image columns per banks.
        maxColID = 69
    elif maxColID <= 35:    #this case for split FOV saphyr with 137 FOVS image columns per bank, but run ID contains 34, 34, 34, or 35 columnIDs (34+34+34+35=137)
        maxColID = 35  

    #n is the number of runIDs per bank.
    n = 0

    #iterate through each runID in runIDList.
    for r in runIDList:
        #if value of runIDMaxColHash[str(r)] equal to maxColID, then reached end of bank. iterate n, add n value to runIDPerBank, and reset n to 0. 
        #runIDPerBank list of number of runIDs for each bank.
        if runIDMaxColHash[str(r)] == maxColID:
            n+=1
            runIDPerBank.append(n)
            n = 0
        
        #if value of runIDMaxColHash[str(r)] less than maxColID, then image columns have not reached end of bank. iterate n and continue.
        elif runIDMaxColHash[str(r)] < maxColID:
            n+=1

    #compute value for runIDPerScan, number of runIDs per scan.
    #There are often issues due to sticking/streaking, low throughput, that result in no dna in some regions of a bank.
    #if this is the case then the value of n (number of runIDs per bank) may be incorrect. This is handled by adding 
    #many differet n values to runIDPerBank. The median value should be correct.
    runIDPerScan = numBanks * np.median(runIDPerBank)     

    return runIDPerScan

#function associates each runID to a scan number.
def groupRunIDByScan(runIDList,runIDPerScan):

    #start with scanNum = 1
    scanNum = 1
    #tempList is a list of runIDs within a single scan number. Will be reset once iterate to next scan.
    tempList = []
    #runIDToScanHash: key is runID value, value is a list with scan number, i.e. [scan#]
    runIDToScanHash = {}

    #iterate through each runID in runIDList.
    for r in runIDList:
        #add runID value to tempList.
        tempList.append(r)

        #if the number of runIDs in tempList equals the value of runIDPerScan (number of runIDs in a scan)
        #then for each runID in tempList, add the scanNum value in a list [scanNum]. Remember, multiple runIDs will likely have the same scan number.
        if len(tempList) == runIDPerScan:
            for t in tempList:
                runIDToScanHash[str(t)] = [scanNum]

            #once the values in tempList have been added to runIDToScanHash, reset tempList to empty, iterate scanNum by 1.
            tempList = []
            scanNum+=1

    #if run is interupted (i.e. did not finish imaging entire flow cell)
    #then len(tempList) will be less runIDPerScan and never be added to runIDToScanHash.
    #if values stilll present in tempList, then fill in remaining values in tempList with scanNum into runIDToScanHash
    if len(tempList) > 0:
        for t in tempList:
            runIDToScanHash[str(t)] = [scanNum]        


    return runIDToScanHash



#function prints how runID and scan numbers were assigned.
def printRunIDInformation(maxColID,runIDPerScan,runIDList,runIDMaxColHash,runIDToScanHash,molNumPerRunID,prefix):
    outfile = open(prefix+'_runID_to_scan.txt','w')
    outfile.write('# maxColID:\t'+str(maxColID)+'\n')
    outfile.write('# Run IDs Per Scan:\t'+str(runIDPerScan)+'\n')

    outfile.write('# RunID\tMaxColumnID\tScanNumber\tNumMoleculesInRunID\n')

    for r in runIDList: #iterate through each r
        outfile.write(str(r)+'\t'+str(runIDMaxColHash[str(r)])+'\t'+str(runIDToScanHash[str(r)][0])+'\t'+str(molNumPerRunID[str(r)])+'\n')

    outfile.close()


#function determines if should split bnx file for each scan number or for give range of scan numbers
#returns scanHash that contains a string of potential Scan## names depending on how bnx file will be parsed.
def parseScanRange(sValue,runIDToScanHash,strScanRange):
    #scanHash. key = scan number, value = Scan## (or if range of scan values selected, Scan##-Scan##)
    scanHash = {}

    #rkeys is list of runID values.
    rkeys = runIDToScanHash.keys()

    if sValue == 1: #if sValue = 1, then generate bnx file for each scan number.
        for r in rkeys:     #iterate through each runID value, r
            s = str(runIDToScanHash[r][0])    #runIDToScanHash[r][0] is scan number

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

#function determines the number of molecules for each bnx file that will be generated
def computeNumMoleculesPerBNX(sValue,runIDList,scanHash,runIDToScanHash,molNumPerRunID):

    #if sValue = 1, then keys will be each scan number, value will be total number of molecules in scan
    #if sValue = 0, then key will be 'one', value will be total number of molecules for given scan range.
    numMolPerFile = {}


    if sValue == 1:     #if sValue = 1, then generate bnx file for each scan number.
        for r in runIDList:
            scanNum = runIDToScanHash[str(r)][0]    #determine scan number value of runID, r

            if str(scanNum) in scanHash:    #if scanNum in scanHash, then add number of molecules for runID, r, to numMolPerFile[str(scanNum)]
                #if only one file being printed (sValue = 1), then numMolPerFile will a key for each scan number
                if str(scanNum) in numMolPerFile:
                    numMolPerFile[str(scanNum)] += molNumPerRunID[str(r)]
                else:
                    numMolPerFile[str(scanNum)] = molNumPerRunID[str(r)]        

    if sValue == 0:     #if sValue = 0, then generate bnx file for range of scan numbers based on user input.
        for r in runIDList:
            scanNum = runIDToScanHash[str(r)][0]    #determine scan number value of runID, r

            if str(scanNum) in scanHash:    #if scanNum in scanHash, then add number of molecules for runID, r, to numMolPerFile['one']
                #if only one file being printed (sValue = 0), then numMolPerFile will have one key, 'one'
                if 'one' in numMolPerFile:
                    numMolPerFile['one'] += molNumPerRunID[str(r)]
                else:
                    numMolPerFile['one'] = molNumPerRunID[str(r)]

    return numMolPerFile




def writeTitle(line,openFileHash,sValue,numMolPerFile):
    skeys = openFileHash.keys()     #list of scan numbers

    if sValue == 0:     #print single bnx file for a range of scan numbers. pick first skey value.

        if line[0:22] == '# Number of Molecules:':  #add the number of molecules to the new bnx file.
            openFileHash[skeys[0]].write('# Number of Molecules:\t'+str(numMolPerFile['one'])+'\n')
        else:
            openFileHash[skeys[0]].write(line)

    elif sValue == 1:     #print bnx file for each scan number.
        for s in skeys:
            if line[0:22] == '# Number of Molecules:':  #add the number of molecules to each scan number bnx file.
                openFileHash[s].write('# Number of Molecules:\t'+str(numMolPerFile[s])+'\n')
            else:
                openFileHash[s].write(line)


#generate hash of new bnx file names.
def generateFileNames(sValue,scanHash,runIDToScanHash,prefix):
    #filenameHash: key=scan number, value is filename.
    filenameHash = {}

    rkeys = runIDToScanHash.keys()

    #if sVavlue=0, then generate single bnx file for range of scan numbers
    #scanHash[s] = ScanRange_1-10_15-20 (for example)
    if sValue == 0:
        for r in rkeys:     #iterate through each runID
            s = str(runIDToScanHash[r][0])    #runIDToScanHash[r][0] is scan number

            if s in scanHash:
                if s in filenameHash:   #if s already in filenameHash, do nothing.
                    donothing = 1
                else:                   #else add s to filenameHash with new file name
                    filenameHash[s] = prefix+'_'+scanHash[s]+'filtered.bnx'

    #if sVavlue=1, then generate bnx file for each scan number
    #scanHash[s] = Scan01, Scan02, etc. for each s.
    elif sValue == 1:
        for r in rkeys:     #iterate through each runID
            s = str(runIDToScanHash[r][0])    #runIDToScanHash[r][0] is scan number

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
def readBNXFile(sValue,numMolPerFile,scanHash,runIDToScanHash,filename,prefix):
    infile = open(filename,'rU')    #open original bnx file. read only!

    #generate hash of new bnx file names based on user input.
    filenameHash = generateFileNames(sValue,scanHash,runIDToScanHash,prefix)
    skeys = filenameHash.keys()     #fkeys is list of 

    #openFiles function: adds open files object to a hash openFileHash
    openFileHash = openFiles(sValue,filenameHash)


    for line in infile: #iterate through each line in bnx file.
        if line[0] == '#':
            writeTitle(line,openFileHash,sValue,numMolPerFile)
        else:
            if line[0] == '0':              #check data from line 0
                keep = 0            #for each molecule, set keep=0, then check runID, scan number to determine if keep=1
                word = line.split() #split line 0 based on white space.

                runID = str(word[11])   #runID, which connects the molecule to the respective scan number
                snum = str(runIDToScanHash[runID][0])    #scan number 

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

#function finds maximum number of imaging columns, runIDs.
maxColID,runIDList,runIDMaxColHash,molNumPerRunID = findNumberColumnsPerBank(args.bnxFile)

#function determines number of runID values per scan.
runIDPerScan = determineRunIDPerScan(maxColID,runIDList,runIDMaxColHash)

#function determines the scan number for each runID.
runIDToScanHash = groupRunIDByScan(runIDList,runIDPerScan)

#print simple text file. Print each scan number, runID value.
printRunIDInformation(maxColID,runIDPerScan,runIDList,runIDMaxColHash,runIDToScanHash,molNumPerRunID,args.prefix)


#Check: If -s 0, code expects a value for -r. If not value give for -r exit code and print warning.
if args.byScan == 0:
    if len(args.scanRange) == 0:
        print
        print
        print '#### WARNING: You have selected to print bnx for range of scan numbers (-s 0), but have not given a range of scan numbers (-r '+args.scanRange+'). This is incompatible. If you select -s 0, then you must give a range of scan numbers (-r). For example: -r 1-10'
        sys.exit()

    elif len(args.scanRange) > 0:   #if values given to -r, continue
        #determine how to parse up the scan numbers based on the range of scan numbers give by user.
        scanHash = parseScanRange(0,runIDToScanHash,args.scanRange)

#Check: If -s 1, then code will print a bnx file for each scan number. There should be no value for -r. 
elif args.byScan == 1:
    if len(args.scanRange) == 0:
        #parse bnx file for each scan number.
        scanHash = parseScanRange(1,runIDToScanHash,args.scanRange)

    elif len(args.scanRange) > 0:   #if values listed for -r, then exit code and print warning.
        print
        print
        print '#### WARNING: You have selected to print a bnx file for each scan number (-s 1), but also given a range of scan numbers (-r '+args.scanRange+'). This is incompatible. If you select -s 1, then you can not select a range of scan numbers (-r) as well.'
        sys.exit()

#compute total number of molecules for each new bnx file.
numMolPerFile = computeNumMoleculesPerBNX(args.byScan,runIDList,scanHash,runIDToScanHash,molNumPerRunID)

print 'Reading, Filtering BNX File:\t',args.bnxFile
#read bnx file and generate new bnx file based on user input.
readBNXFile(args.byScan,numMolPerFile,scanHash,runIDToScanHash,args.bnxFile,args.prefix)











