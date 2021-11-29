from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from optparse import OptionParser

import numpy as np
import pandas as pd
import seaborn as sns
import os, multiprocessing, progressbar, signal, time


def init_worker(tqdm_lock=None):
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    if tqdm_lock is not None:
        tqdm.set_lock(tqdm_lock)

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
    cds_validity.to_csv(f'./data/table/{org}_{subset}_seq_validity.csv', sep=',', index=False)
    return cds_validity


def get_nr_mutants(cds_id, cds):
    set_mutations = {
                        'A': ['T', 'C', 'G'],
                        'T': ['A', 'C', 'G'],
                        'C': ['A', 'G', 'T'],
                        'G': ['A', 'C', 'T']
                    }
    all_mut_dict = {}
    all_mut_dict[f'{cds_id}_WT'] = str(cds.translate()).replace('*','')
    for i in range(len(cds)) :
        for j in range(3):
            mutant = cds.tomutable()
            mutant[i] = set_mutations[cds[i]][j]
            seq = str(mutant.toseq().translate()).replace('*', '')
            if seq not in all_mut_dict.values():
                all_mut_dict[f'{cds_id}_{i}_{cds[i]}{set_mutations[cds[i]][j]}'] = str(mutant.toseq().translate()).replace('*', '')
    return all_mut_dict


def createSeqFile(args):
    seq_id = args[0]
    seqDict = args[1]
    outPath = args[2]

    cds_id = seq_id
    cds = seqDict[seq_id].seq
    folder = f'{outPath}/{cds_id}'
    if not os.path.isdir(folder) :
        os.mkdir(folder)

    f = open(f'{folder}/{cds_id}.seq', 'w')
    mutation_dict = get_nr_mutants(cds_id, cds)
    for key in mutation_dict.keys():
        f.write(f'{key} N N 7 298 0.1 {mutation_dict[key]}\n')
    f.close()


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-o", "--organism", dest="org", default="None", help="[Required] Provide the name of the organism (MM or HG)")
    parser.add_option("-s", "--subset", dest="subset", default="None", help="[Required] Provide subset name (chap_client)")
    parser.add_option("-e", "--output_path", dest="outPath", default="None", help="[Required] Provide name of the result folder")
    (options, args) = parser.parse_args()
    org = options.org
    subset = options.subset
    outPath = options.outPath

    myfolders =  [ os.path.join(outPath, folder) for folder in os.listdir(outPath) ]
    all_agg_scores = pd.read_csv('../../data/aggregation_propensity/HGMM_agg_scores.csv')
    print(f'Number of proteins to mutate: {len(all_agg_scores)}\n')

    #### Collecting cds for these proteins
    print(f'Collecting CDS for {org}\n')
    chap_clt_sub = chap_clt_subset(all_agg_scores)
    if 'MM' in org :
        fastaFile = '../../data/ortholog_dataset/uni_MM_cds_orthologs.faa'
        if 'chap_client' in subset:
            idList = chap_clt_sub['proteinID_y'].values
        else :
            subset ='others'
            idList = all_agg_scores[~all_agg_scores['proteinID_y'].isin(chap_clt_sub['proteinID_y'])].values
    if 'HG' in org:
        fastaFile = '../../data/ortholog_dataset/uni_HG_cds_orthologs.faa'
        if 'chap_client' in subset:
            idList = chap_clt_sub['proteinID_x'].values
        else :
            subset = 'others'
            idList = all_agg_scores[~all_agg_scores['proteinID_x'].isin(chap_clt_sub['proteinID_x'])].values

    seqSubset = collect_subset(fastaFile, idList)
    print(f'Number of proteins collected: {len(seqSubset)}')

    #### Checking cds_validity
    cds_validity = isValid(seqSubset, org, subset)
    valid_seq = cds_validity[cds_validity['valid_cds'] == True]['proteinID'].values
    print(f'Number of valid proteins: {len(valid_seq)}\n')

    #### Create seq file
    print('Creation of Seq files')
    args =  [ (seq_id, seqSubset, outPath) for seq_id in valid_seq ]
    p = multiprocessing.Pool(initializer=init_worker, initargs=(tqdm.get_lock(),), processes=10)
    try:
        pbar = tqdm(myfolders, maxinterval=1.0, miniters=1, desc="Written files", bar_format="{desc}:{percentage:3.0f}%|{bar}|")
        for _, result in enumerate(p.imap_unordered(createSeqFile, args, chunksize=1)):
            # if result:
            #     print("File already written !")
            pbar.update(1)  # Everytime the iteration finishes, update the global progress bar
        pbar.close()
        p.close()
        p.join()
    except KeyboardInterrupt:
        print("KeyboardInterrupt, terminating workers.")
        pbar.close()
        p.terminate()
        p.join()
        exit(1)
