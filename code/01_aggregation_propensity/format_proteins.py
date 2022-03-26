#!/usr/bin/env python

'''
Author: Savandara Besse
Institution: University of Montreal
Created: 12-21-2018

Change U aa to X
'''
import progressbar, re

from Bio import SeqIO
from Bio.Seq import Seq
from optparse import OptionParser


def modifySequences(fastaFile):
    newRecords = []
    for seqRecord in SeqIO.parse(fastaFile, format='fasta'):
        seqRecord.name = ''
        seqRecord.description = ''
        if 'U' in seqRecord.seq :
            seqRecord.seq = Seq(str(seqRecord.seq).replace('U','X'))
            newRecords.append(seqRecord)
        else:
            newRecords.append(seqRecord)
    return newRecords

def saveSequences(outputFile, sequences):
    f = open(outputFile, 'w')
    SeqIO.write(sequences, f,'fasta')
    f.close()

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-i", "--input_fasta", dest="fastaFile", default="None", help="[Required] Provide the path of the initial FASTA file")
    parser.add_option("-o", "--output_fasta", dest="outputFile", default="None", help="[Required] Provide the FASTA destination file ")
    (options, args) = parser.parse_args()
    fastaFile = options.fastaFile
    outputFile = options.outputFile

    sequences = modifySequences(fastaFile)
    saveSequences(outputFile, sequences)
