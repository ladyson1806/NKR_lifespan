#!/usr/bin/env python

'''
Author: Savandara Besse
Created: 02-14-2018
Modified: 04-11-2019

Parse sql table output from Inparalog software to create
separated fasta files of ortholog proteins for each studied organism

Modification:
Change way of recording FASTA sequences
'''

import pandas as pd
import os, progressbar, re, threading
from Bio import SeqIO
from optparse import OptionParser


def splitTable(dataframe, columnName, columnValue):
    return dataframe[dataframe[columnName] == columnValue]

def getTable(table, header, columnName, columnValue):
    df = pd.read_table(table, header=None)
    df.columns = header
    df2 = splitTable(df, columnName, columnValue)
    return df2

def getRecord(initFasta, ID):
    for seqRecord in SeqIO.parse(initFasta, "fasta"):
        if '|' in seqRecord.id :
            toSearch = seqRecord.id.split('|')[1]
        else:
            toSearch = seqRecord.id
        try:
            if ID == toSearch:
                return seqRecord
        except:
            print(ID)
            return None

def saveProtein_blasdtdbcmd(proteome, ID, orthoFile):
    os.system('blastdbcmd -db {} -entry {} >> {}'.format(proteome, ID, orthoFile))

def createFastaFile(dataframe, columnName, proteome, outputFile):
    bar = progressbar.ProgressBar()
    for ID in bar(dataframe[columnName].values):
        saveProtein_blasdtdbcmd(proteome, ID, outputFile)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-x", "--organism", dest="organism", default="None", help="[Required] Organism name as it is written in sql table")
    parser.add_option("-p", "--proteome", dest="fasta", default="None", help="[Required] Protein fasta file")
    parser.add_option("-t", "--table", dest="table", default="None", help="[Required] Ortholog sql table from Inparanoid")
    parser.add_option("-o", "--destinationFile", dest="output", default="None", help="[Required] Name of the output table")
    (options, args) = parser.parse_args()
    organism = options.organism
    initTable = options.table
    proteome = options.fasta
    output = options.output

    header = ['clusterNumber', 'bestScore', 'Organism', 'Score', 'proteinID', 'Bootstrap']
    subTable = getTable(initTable, header,'Organism', organism)
    createFastaFile(subTable, 'proteinID', proteome, output)
