#!/usr/bin/env python

import sys
import os
import csv
import subprocess
from shutil import copy2
from collections import namedtuple

from pbcore.io.FastaIO import FastaReader, FastaWriter, FastaRecord

BlasrM4 = namedtuple('BlasrM4', ['qname', 'tname', 'score', 'pctsimilarity', 
                                 'qstrand', 'qstart', 'qend', 'qseqlength',
                                 'tstrand', 'tstart', 'tend', 'tseqlength',
                                 'mapqv'])

MIN_LENGTH = 500
MIN_SIMILARITY = 0.9
NUM_PROC = 1

def fileExists( filename ):
    return os.path.exists(filename) and (os.path.getsize(filename) > 0)

def isExe( filePath ):
    if filePath is None:
        return False
    return os.path.isfile(filePath) and os.access(filePath, os.X_OK)

def which(program):
    """
    Find and return path to local executables  
    """
    fpath, fname = os.path.split(program)
    if fpath:
        if isExe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exeFile = os.path.join(path, program)
            if isExe(exeFile):
                return exeFile
    return None

def validateInputFile( fileName, allowedSuffixes ):
    allowedSuffixes = [allowedSuffixes] if isinstance(str, type(allowedSuffixes)) \
                                        else allowedSuffixes
    # First we check whether the input file has a valid suffix
    try:            
        assert any( [fileName.endswith(suffix) for suffix in allowedSuffixes] )
    except AssertionError:
        raise ValueError('File does not have an allowed suffix! %s' % \
                                                        allowedSuffixes)
    # Next we check whether the input file exists where specified
    try: 
        assert fileExists( fileName )
    except AssertionError:
        raise OSError('Input file does not exist!')
    # Finally we return the absolute path to the file
    return os.path.abspath( fileName )

def validateExecutable( executable ):
    exePath = which( executable )
    try: 
        assert exePath is not None
    except AssertionError:
        raise ValueError('"%s" is not a valid executable!' % executable)
    return exePath

def validateInt( integer, minValue=None, maxValue=None ):
    try: # First we check whether the supplied parameter is an Int
        assert isinstance(integer, int)
    except AssertionError:
        raise TypeError('Parameter is not an Integer!')
    # If a minimum value is supplied, compare it to the Integer
    if minValue is not None:
        try:
            assert integer >= minValue
        except AssertionError:
            raise ValueError("Integer is less than Minimum Value!")
    # If a maximum value is supplied, compare it to the Integer
    if maxValue is not None:
        try:
            assert integer <= maxValue
        except AssertionError:
            raise ValueError("Integer is greater than Maximum Value!")
    return integer

def validateFloat( floating_point, minValue=None, maxValue=None ):
    try: # First we check whether the supplied parameter is an Int
        assert isinstance(floating_point, float)
    except AssertionError:
        raise TypeError('Parameter is not a Floating Point!')
    # If a minimum value is supplied, compare it to the Integer
    if minValue is not None:
        try:
            assert floating_point >= minValue
        except AssertionError:
            raise ValueError("Float is less than Minimum Value!")
    # If a maximum value is supplied, compare it to the Integer
    if maxValue is not None:
        try:
            assert floating_point <= maxValue
        except AssertionError:
            raise ValueError("Float is greater than Maximum Value!")
    return floating_point

class ContaminationRemover(object):
    """
    A tool to remove vector and E. coli contamination from 
    assembled or pre-assembled reads
    """

    def __init__(self, sequenceFile=None, vectorFile=None, genomeFile=None):
        if sequenceFile is None:
            self.initializeFromArgs()
        else:
            self.initializeFromCall( sequenceFile, vectorFile, genomeFile )
        self.validateSettings()

    def initializeFromArgs(self):
        import argparse
        desc = 'A tool to remove vector and E. coli contamination from' + \
                ' assembled or pre-assembled reads'
        parser = argparse.ArgumentParser(description=desc)
        parser.add_argument('sequenceFile', metavar='SEQUENCE_FILE',
                            help='File of sequence data to filter')
        parser.add_argument('--vector', metavar='VECTORS',
                            dest='vectorFile',
                            help='File of vector sequences to remove')
        parser.add_argument('--genome', metavar='GENOME',
                            dest='genomeFile',
                            help='File of bacterial genomes to remove')
        parser.add_argument('--num_proc', default=NUM_PROC,
                            dest='numProc', type=int, metavar='INT',
                            help='Number of processors to use')
        parser.add_argument('--min_similarity', default=MIN_SIMILARITY,
                            dest='minSimilarity', type=float, metavar='FLOAT',
                            help='Minimum percent similarity for Blasr hits')
        parser.add_argument('--min_length', default=MIN_LENGTH,
                            dest='minLength', type=int, metavar='INT',
                            help='Minimum length to allow in filtered reads')
        parser.add_argument('--blasr', metavar='EXE_PATH', default='blasr',
                            help='Blasr executable to use for comparisons')
        args = parser.parse_args()
        self.__dict__.update( vars(args) )

    def initializeFromCall( self, sequenceFile, vectorFile, genomeFile ):
        self.sequenceFile = sequenceFile
        self.vectorFile = vectorFile
        self.genomeFile = genomeFile

    def validateSettings( self ):
        fastaSuffix = ['.fa', '.fasta']
        # Validate files
        self.sequenceFile = validateInputFile( self.sequenceFile, fastaSuffix )
        if self.vectorFile is not None:
            self.vectorFile = validateInputFile( self.vectorFile, fastaSuffix )
        if self.genomeFile is not None:
            self.genomeFile = validateInputFile( self.genomeFile, fastaSuffix )
        # Validate Executables
        self.blasr = validateExecutable( self.blasr )
        # Validate other variables
        validateInt( self.numProc, minValue=1 )
        validateInt( self.minLength, minValue=0 )
        validateFloat( self.minSimilarity, minValue=0.0, maxValue=1.0)
        # Finally initialize a process counter and bin for temp files
        self.counter = 0
        self.tempFiles = []

    def runBlasr( self, queryFile, referenceFile ):
        print "Beginnging Blasr analysis Round #%s..." % self.counter
        print "Running Blasr..."
        outputFile = 'temp_%s.out' % self.counter
        p = subprocess.Popen( [self.blasr,
                               '-bestn', '1',
                               '-nproc', str(self.numProc),
                               '-m', '4',
                               '-out', outputFile,
                               queryFile,
                               referenceFile] )
        p.wait()
        self.tempFiles.append( outputFile )
        return outputFile

    def parseBlasrOutput( self, blasrOutput ):
        print "Parsing the Blasr output..."
        blasrHits = {}
        with open( blasrOutput, 'r') as handle:
            for row in csv.reader( handle, delimiter=' ' ):
                hit = BlasrM4._make( row )
                if hit.pctsimilarity < self.minSimilarity:
                    continue
                rec_id = hit.qname.split('/')[0]
                blasrHits[rec_id] = hit
        return blasrHits

    def readSequenceData( self, sequenceFile ):
        print "Reading sequence data into memory..."
        sequenceData = {}
        for record in FastaReader( sequenceFile ):
            sequenceData[ record.name ] = record
        return sequenceData

    def trimSequenceData( self, sequenceData, blasrHits ):
        print "Trimming out vector sequence..."
        trimmedSeqData = []
        for rec_id, record in sequenceData.iteritems():
            # If the record has a Blasr hit, find it
            try:
                hit = blasrHits[rec_id]
            # Otherwise keep the sequence as-is
            except KeyError:
                trimmedSeqData.append(record)
                continue
            # For records with hits, cut out and keep the good sequence
            start = int(hit.qstart)
            end = int(hit.qend)
            if start > self.minLength:
                newName = record.name + '_5p'
                newSequence = record.sequence[:start]
                newRecord = FastaRecord( newName, newSequence )
                trimmedSeqData.append( newRecord )
            if len(record.sequence) - end > self.minLength:
                newName = record.name + '_3p'
                newSequence = record.sequence[end:]
                newRecord = FastaRecord( newName, newSequence )
                trimmedSeqData.append( newRecord )
        return trimmedSeqData

    def removeContaminatedSequences( self, sequenceData, blasrHits ):
        print "Trimming out vector sequence..."
        filteredSeqData = []
        for rec_id, record in sequenceData.iteritems():
            # If the record has a Blasr hit, find it
            try:
                hit = blasrHits[rec_id]
            # If not, we keep the record
            except KeyError:
                filteredSeqData.append(record)
        return filteredSeqData

    def writeSequenceData( self, sequenceData ):
        outputFile = 'temp_%s.fasta' % self.counter
        with FastaWriter( outputFile ) as handle:
            for record in sequenceData:
                handle.writeRecord( record )
        self.tempFiles.append( outputFile )
        return outputFile

    def removeTempFiles( self ):
        print "Cleaning up temporary files.."
        for filename in self.tempFiles:
            os.remove( filename )

    def filterGenomicSequence( self, inputFile ):
        print "Filtering out reads with Genomic contamination"
        self.counter += 1
        blasrOutput = self.runBlasr( inputFile, 
                                     self.genomeFile )
        blasrHits = self.parseBlasrOutput( blasrOutput )
        if not blasrHits:
            return inputFile
        sequenceData = self.readSequenceData( inputFile )
        filteredSequences = self.removeContaminatedSequences( sequenceData,
                                                              blasrHits )
        filteredFile = self.writeSequenceData( filteredSequences )
        return filteredFile

    def filterVectorSequence( self, sequenceFile ):
        print "Trimming vector sequence from sequence reads"
        currSeqFile = sequenceFile
        while True:
            self.counter += 1
            blasrOutput = self.runBlasr( currSeqFile, 
                                         self.vectorFile )
            blasrHits = self.parseBlasrOutput( blasrOutput )
            # If there are no good Blasr Hits
            if not blasrHits:
                break
            # Trim the identified regions from
            sequenceData = self.readSequenceData( currSeqFile )
            trimmedData = self.trimSequenceData( sequenceData, blasrHits )
            currSeqFile = self.writeSequenceData( trimmedData )
        return currSeqFile

    def backupInputFile( self ):
        path, filename = os.path.split( self.sequenceFile )
        print 'Backing up "%s"...' % filename
        backupFile = 'BACKUP_' + filename
        backupPath = os.path.join(path, backupFile)
        if os.path.exists( backupPath ):
            raise IOError('Backup of "%s" already exists!' % filename )
        copy2( self.sequenceFile , backupPath )

    def renameOutputFile( self, outputFile ):
        print "Renaming output file..."
        copy2( outputFile, self.sequenceFile )

    def __call__(self):
        sequenceFile = self.sequenceFile
        # If a genome file was specified, remove any contaminated reads
        if self.genomeFile:
            sequenceFile = self.filterGenomicSequence( sequenceFile )
        # If a vector file was specified, trim away any vector sequence
        if self.vectorFile:
            sequenceFile = self.filterVectorSequence( sequenceFile )
        # If no changes
        if sequenceFile == self.sequenceFile:
            print "No filtering or trimming required, skipping backup..."
        else:
            self.backupInputFile()
            self.renameOutputFile( sequenceFile )
        # Clean-up
        self.removeTempFiles()
        print "Process Complete"

if __name__ == '__main__':
    remover = ContaminationRemover()
    remover()
