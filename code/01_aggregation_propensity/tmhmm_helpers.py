import re
import pandas as pd

keywords = {
    'length': 'Length',
    'TMHs': 'Number of predicted TMHs',
    'TMHs_aa': 'Exp number of AAs in TMHs',
    '60st_aa': 'Exp number, first 60 AAs',
    'N_prob': 'Total prob of N-in',
}


def returnKey(sentence):
    for key in keywords.keys():
        if sentence in keywords[key]:
            return key


def collectValue(ID, line, Dict):
    if ID not in Dict.keys():
        Dict[ID] = {}
        Dict[ID]['Prob_N_sign'] = 'No'
    if 'POSSIBLE N-term signal sequence' in line:
        Dict[ID]['Prob_N_sign'] = 'Yes'
    else:
        if 'Length' in line:
            key = returnKey('Length')
        elif 'Number of predicted TMHs' in line:
            key = returnKey('Number of predicted TMHs')
        elif 'Exp number of AAs in TMHs' in line:
            key = returnKey('Exp number of AAs in TMHs')
        elif 'Exp number, first 60 AAs' in line:
            key = returnKey('Exp number, first 60 AAs')
        elif 'Total prob of N-in' in line:
            key = returnKey('Total prob of N-in')
        value = re.findall(r'([0-9.]+)$', line)[0]
        Dict[ID][key] = float(value)
    return Dict


def collect_TMHMM_results(TMHMM_file):
    f = open(TMHMM_file,'r')
    refComments = [ x.strip() for x in f.readlines()]
    f.close()

    tmp = {}
    for comment in refComments :
        if '#' in comment:
            ID = comment.split(' ')[1]
            ref = collectValue(ID, comment, tmp)

    TM_table = pd.DataFrame.from_dict(tmp, orient='index')
    TM_table = TM_table.reset_index()
    TM_table = TM_table.rename(index=str, columns={"index": "proteinID"})
    TM_table = TM_table[['proteinID','length', 'TMHs_aa', '60st_aa', 'N_prob', 'TMHs', 'Prob_N_sign']]
    TM_table['TM_protein'] = TM_table['TMHs'] >= 1
    return TM_table
