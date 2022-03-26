import pandas as pd
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
import json, os, progressbar, re, time

from Bio import SeqIO
from Bio import AlignIO

def return_pfam_entry(proteinID):
    for i in range(len(data['results'])):
        if data['results'][i]['metadata']['accession'] == proteinID :
            return data['results'][i]
        
def get_domain_nb(proteinID):
    json_result = return_pfam_entry(proteinID)
    cpt = 0
    for pfam_entry in json_result['entry_subset']:
        cpt += len(pfam_entry['entry_protein_locations'])
    return cpt

def collect_nb_dom(MM_IDs):
    tmp = []
    bar = progressbar.ProgressBar()
    for ID in bar(MM_IDs):
        try:
            tmp.append((ID, get_domain_nb(ID)))
        except:
            tmp.append((ID, 0))

    nb_dom = pd.DataFrame(tmp, columns=['proteinID_y', 'nb_domains'])
    return nb_dom

def get_sequences(y):
    tmp = []
    
    x = ortho_pairs[ortho_pairs['proteinID_y'] == y]['proteinID_x'].values[0]
    for seqRecord in SeqIO.parse(MM_fasta, format='fasta'):
        if y in seqRecord.id : 
            tmp.append(seqRecord)
    for seqRecord in SeqIO.parse(HG_fasta, format='fasta'):
        if x in seqRecord.id :
            tmp.append(seqRecord)
    return tmp

def build_alignment(y):
    records = get_sequences(y)
    with open("/media/savvy/DATA3/savvy/project_2018/MSA/seq_tmp.fasta", "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
    if os.path.exists(f'./tmp/{records[0].id}_{records[1].id}.fasta'):
        pass 
    else:
        os.system(f'~/bin/muscle -in /media/savvy/DATA3/savvy/project_2018/MSA/seq_tmp.fasta -out /media/savvy/DATA3/savvy/project_2018/MSA/{records[0].id}_{records[1].id}.fasta -maxiters 1 -diags -sv -distance1 kbit20_3')


HG_fasta = './data/ortholog_dataset/uni_HG_orthologs.faa'
MM_fasta = './data/ortholog_dataset/uni_MM_orthologs.faa'

MM_IDs = [ seqRecord.id for seqRecord in SeqIO.parse(MM_fasta, format='fasta')]

ortho_pairs = pd.read_csv('./data/TABLE/HGMM_agg_scores.csv', sep=',') 


bar = progressbar.ProgressBar()
for ID in bar(ortho_pairs['proteinID_y']):
    build_alignment(ID)
