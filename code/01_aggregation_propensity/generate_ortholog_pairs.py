#!/usr/bin/env python

'''
Author: Savandara Besse
Created: 03-15-2018

Collect the best pair of orthologs
for each inparanoid cluster and build a dataframe
to save the candidate proteins
'''

import progressbar
import pandas as pd

def splitTable(dataframe, columnName, columnValue):
    return dataframe[dataframe[columnName] == columnValue]

def getTable(table, header, columnName, columnValue):
    df = pd.read_table(table, header=None)
    df.columns = header
    df2 = splitTable(df,columnName, columnValue)
    return df2

def collectOrthologs(df):
    orthoList = []
    bar = progressbar.ProgressBar()
    nbCluster = df['clusterNumber'].unique().tolist()
    for cluster in bar(nbCluster):
        tmp = df[df['clusterNumber'] == cluster].values.tolist()[0]
        if tmp not in orthoList:
            orthoList.append(tmp)
    return orthoList

def buildOrthologTable(df, orthoList):
    columns = df.columns.tolist()
    dff = pd.DataFrame(data=orthoList, columns=columns)
    dff.to_csv('./uniprot_results/HS_MM_uniOrthologGroups.csv',sep='\t', index=False)

if __name__ == '__main__':
    # sqltable='/media/savvy/DATA3/savvy/inparanoid/HG-MM_Uniprot/sqltable.MM-HG'
    # header = ['clusterNumber', 'bestScore', 'Organism', 'Score', 'proteinID', 'Bootstrap']
    # HGTable = getTable(sqltable, header,'Organism', 'HG')
    # MMTable = getTable(sqltable, header,'Organism', 'MM')

    sqltable='/media/savvy/DATA3/savvy/inparanoid/HS-MM_Inparanoid/sqltable.H.sapiens-M.musculus'
    header = ['clusterNumber', 'bestScore', 'Organism', 'Score', 'proteinID', 'Bootstrap']
    HGTable = getTable(sqltable, header,'Organism', 'H.sapiens')
    MMTable = getTable(sqltable, header,'Organism', 'M.musculus')



    ## Collect pair-wise orthologs
    df = HGTable.merge(MMTable, on='clusterNumber')
    ## Collect unique ortholog groups
    orthoList = collectOrthologs(df)
    buildOrthologTable(df, orthoList)
