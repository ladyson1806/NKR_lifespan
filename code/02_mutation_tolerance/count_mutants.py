from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Seq import MutableSeq
from optparse import OptionParser

import numpy as np
import pandas as pd
import seaborn as sns
import progressbar, time

def chap_clt_subset(all_agg_scores):
    uniprot_mapping = pd.read_csv('../../data/chaperone_clients/human_ensembl_to_uniprot.tab', sep='\t')
    hs_mm_orthologs = pd.read_csv('../../data/chaperone_clients/HS_MM_uni_ortholog_groups.csv', sep='\t')
    hs_mm_orthologs = hs_mm_orthologs[['proteinID_x', 'proteinID_y']]
    mm_chap_clt = hs_mm_orthologs[hs_mm_orthologs['proteinID_x'].isin(uniprot_mapping['Entry'])]['proteinID_y']
    chap_clt_sub = all_agg_scores[all_agg_scores['proteinID_y'].isin(list(mm_chap_clt))]
    return chap_clt_sub


def collect_subset(fastaFile, idList):
    seqDict = {}
    for seqRecord in SeqIO.parse(fastaFile, 'fasta'):
        if seqRecord.id in idList :
            seqDict[seqRecord.id] = seqRecord
    return seqDict


def isValid(seqDict, org, subset):
    tmp = []
    bar = progressbar.ProgressBar()
    for key in bar(seqDict.keys()):
        cds_id = key
        cds = seqDict[key].seq
        if len(cds) % 3 == 0 :
            if list(np.unique(cds)) != ['A', 'C', 'G', 'T']:
                tmp.append([cds_id, 'non_standard_nucleotide', False])
            else:
                tmp.append([cds_id, 'standard', True])
        else:
            tmp.append([cds_id, 'non_standard_length', False])
    cds_validity = pd.DataFrame(tmp, columns=['proteinID', 'description', 'valid_cds'])
    # cds_validity.to_csv(f'../../data/mutation_tolerance/{org}_{subset}_seq_validity.csv', sep=',', index=False)
    return cds_validity


def mutate_seq(cds_id, cds):
    set_mutations = {
                        'A': ['T', 'C', 'G'],
                        'T': ['A', 'C', 'G'],
                        'C': ['A', 'G', 'T'],
                        'G': ['A', 'C', 'T']
                    }

    all_mut_dict = {}
    wt = str(cds.translate(to_stop=True))
    all_mut_dict[f'{cds_id}_WT'] = str(cds.translate(to_stop=True))
    for i in range(len(cds)) :
        for j in range(3):
            mutant = MutableSeq(str(cds))
            mutant[i] = set_mutations[cds[i]][j]
            mt = str(Seq(str(mutant)).translate(to_stop=True))
            if (len(wt) != len(mt)) or (mt[0] != 'M'):
                pass 
            else :
                all_mut_dict[f'{cds_id}_{i}_{cds[i]}{set_mutations[cds[i]][j]}'] = mt
    return all_mut_dict


def mutant_table(valid_seq, seqDict, org, subset):
    bar = progressbar.ProgressBar()
    all_seq_info = []
    for seq_id in bar(valid_seq) :
        cds_id = seq_id
        cds = seqDict[seq_id].seq
        mutation_dict = mutate_seq(cds_id, cds)
        all_seq_info.append([cds_id, len(mutation_dict.values()), len(set(mutation_dict.values()))])
    compMut_info = pd.DataFrame(all_seq_info, columns=['proteinID', 'all_MT', 'unique_MT'])
    compMut_info.to_csv(f'../../data/mutation_tolerance/{org}_{subset}_mutants_counts.csv', sep=',', index=False)
    print(compMut_info)


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-o", "--organism", dest="org", default="None", help="[Required] Provide the name of the organism (MM or HG)")
    parser.add_option("-s", "--subset", dest="subset", default="None", help="[Required] Provide subset name (all or atx or chap_client)")
    (options, args) = parser.parse_args()
    org = options.org
    subset = options.subset

    all_agg_scores = pd.read_csv('../../data/aggregation_propensity/HGMM_agg_scores.csv')
    print(f'Number of proteins to mutate: {len(all_agg_scores)}\n')

    #### Collecting cds for these proteins
    print(f'Collecting CDS for {org}\n')
    if 'MM' in org :
        fastaFile = '../../data/ortholog_dataset/uni_MM_cds_orthologs.faa'
        if 'chap_client' in subset:
            chap_clt_sub = chap_clt_subset(all_agg_scores)
            idList = chap_clt_sub['proteinID_y'].values
        elif 'atx' in subset:
            idList = ['P28658', 'Q9CVD2']
        else :
            subset ='all'
            idList = all_agg_scores['proteinID_y'].values
    if 'HG' in org:
        fastaFile = '../../data/ortholog_dataset/uni_HG_cds_orthologs.faa'
        if 'chap_client' in subset:
            chap_clt_sub = chap_clt_subset(all_agg_scores)
            idList = chap_clt_sub['proteinID_x'].values
        elif 'atx' in subset:
            idList = ['G5BVC0', 'G5AZL7']
        else :
            subset = 'all'
            idList = all_agg_scores['proteinID_x'].values

    seqSubset = collect_subset(fastaFile, idList)
    print(f'Number of proteins collected: {len(seqSubset)}')

    #### Checking cds_validity
    cds_validity = isValid(seqSubset, org, subset)
    valid_seq = cds_validity[cds_validity['valid_cds'] == True]['proteinID'].values
    print(f'Number of valid proteins: {len(valid_seq)}\n')

    # #### Create mutant table for all or chap_client subset
    # print('Creation of mutant table')
    # mutant_table(valid_seq, seqSubset, org, subset)

    #### ATX proteins
    mutant_table(valid_seq, seqSubset, org, subset)
