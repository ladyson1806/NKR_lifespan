#!/usr/bin/env python

'''
Author: Savandara BESSE
Institution: University of Montreal
Created: 03-14-2018
Modified: 04-11-2019
'''

import os
import pandas as pd
import in_getOrthologs as getOrth
from Bio import SeqIO

def buildTableWithLength(df, idColumn, initFasta, lengthColumn):
    orgTable = df[['clusterNumber', idColumn]]
    length = []
    for seqRecord in SeqIO.parse(initFasta,'fasta'):
        if '|' in seqRecord.id :
            ID = seqRecord.id.split('|')[1]
        else :
            ID = seqRecord.id
        if ID in orgTable[idColumn].values:
            length.append([ID, len(seqRecord.seq)])
    tmp = pd.DataFrame(data=length, columns=[idColumn, lengthColumn])
    lengthTable = orgTable.merge(tmp, on=idColumn)
    return lengthTable

if __name__ == '__main__':
    print('Collecting best ortholog pairs')
    os.system('python3 groupOrthologs.py')
    df = pd.read_csv('./uniprot_results/HG_MM_uniOrthologGroups.csv', sep='\t')

    ## Collect onrthologs with a good bootstrapping Score
    print('Collecting best ortholog pairs with high bootstrapping score')
    df = df[(df['Bootstrap_x'] == '100%') & (df['Bootstrap_y'] == '100%')]
    print('Final number of ortholog pairs: ', len(df))

    ## Collect length for HG and MM proteins
    HG_Orth = buildTableWithLength(df, 'proteinID_x', './uniprot_results/HG_orthologs.faa', 'length_x')
    MM_Orth = buildTableWithLength(df, 'proteinID_y', './uniprot_results/MM_orthologs.faa', 'length_y')

    ## Obtain final dataframe with HG and MM protein length
    dff = HG_Orth.merge(MM_Orth, on='clusterNumber')
    dff.to_csv('./uniprot_results/HG_MM_Orthologs_Length.csv', sep='\t', index=False)

    # ## Create fasta files
    getOrth.createFastaFile(dff, 'proteinID_x', '/media/DATA1/savvy/BLAST/model/het_glab/het_glab', './uniprot_results/uni_HG_orthologs.faa' )
    getOrth.createFastaFile(dff, 'proteinID_y', '/media/DATA1/savvy/BLAST/model/mus_musc/mus_musc', './uniprot_results/uni_MM_orthologs.faa')

    # # ## Create reduced fasta files
    # dff = dff[ (dff['length_x'] < 400 ) & (dff['length_y'] < 400) & (dff['length_x'] > 300 ) & (dff['length_y'] > 300)] ## Range decided based on median length for eukaryote organisms
    # dff.to_csv('./custom_results/reduced_HG_MM_Orthologs.csv', sep='\t', index=False)
    # getOrth.createFastaFile('./custom_results/reduced_HG_orthologs.faa', dff, 'proteinID_x', HG_AllOrth)
    # getOrth.createFastaFile('./custom_results/reduced_MM_orthologs.faa', dff, 'proteinID_y', MM_AllOrth)
