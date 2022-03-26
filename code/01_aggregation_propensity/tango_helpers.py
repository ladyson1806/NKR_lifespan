import pandas as pd
import os, progressbar
from Bio import SeqIO


def collectSeq(proteome):
    seqCollection = []
    for seqRecord in SeqIO.parse(proteome, format='fasta'):
        seqCollection.append((seqRecord.id, str(seqRecord.seq)))
    return seqCollection

def createSEQ(resPath, seqTuple):
    f = open(os.path.join(resPath,'{}.seq'.format(seqTuple[0])), mode='w')
    f.write('{} N N 7 298 0.1 {}\n'.format(seqTuple[0], seqTuple[1]))
    f.close()

def collect_aggScore(tangoTable):
    proteinID = tangoTable.split('/')[-1].replace('.txt', '')
    tmp = pd.read_csv(tangoTable, sep='\t')
    agg_score = sum(tmp['Aggregation']) / len(tmp)
    return proteinID, agg_score

def collect_aggTable(tangoTable):
    proteinID = tangoTable.split('/')[-1].replace('.txt', '')
    tmp = pd.read_csv(tangoTable, sep='\t')
    return tmp