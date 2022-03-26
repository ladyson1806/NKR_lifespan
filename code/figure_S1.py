import json, os, re

import numpy as np
import pandas as pd
import plotnine as p9
import scipy.stats as stats

from Bio import SeqIO
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests
from figure_3 import generate_input_for_figure

def check_annotations(ID, fastaFile):
    for seqRecord in SeqIO.parse(fastaFile, format='fasta'):
        if ID in seqRecord.id :
            if '_MOUSE' in seqRecord.description:
                description = re.findall(r'_MOUSE (.*) OS', seqRecord.description)[0]
            if '_HETGA' in seqRecord.description:
                description = re.findall(r'_HETGA (.*) OS', seqRecord.description)[0]
            return description

def chaperone_clients_subset():
    uniprot_mapping = pd.read_csv('../data/chaperone_clients/human_ensembl_to_uniprot.tab', sep='\t')
    hs_mm_orthologs = pd.read_csv('../data/chaperone_clients/HS_MM_uni_ortholog_groups.csv', sep='\t')
    hs_mm_orthologs = hs_mm_orthologs[['proteinID_x', 'proteinID_y']]
    mm_chap_clt = hs_mm_orthologs[hs_mm_orthologs['proteinID_x'].isin(uniprot_mapping['Entry'])]['proteinID_y']
    return mm_chap_clt


def build_gene_list(mm_chap_clt):                                                                                                                                                                      
    MM_fasta = '../data/ortholog_dataset/uni_MM_orthologs.faa'

    main_path = os.path.dirname(os.getcwd())     
    dom_path = os.path.join(main_path, 'data/over-representation_analysis/domains/OUTPUTS/JSON/')
    seq_path = os.path.join(main_path, 'data/over-representation_analysis/whole_sequence/OUTPUTS/JSON/')

    all_gene_list = []
    for json_path in [seq_path, dom_path]:
        if 'domains' in json_path:
            analysis = 'Domains'
        elif 'whole_sequence' in json_path :
            analysis = 'Proteins'
    
        for file in os.listdir(json_path):
            with open(os.path.join(json_path, file), 'r') as j:
                contents = json.loads(j.read())
                if 'high' in file :
                    phenotype = 'Higher aggregation propensity in NKM'
                elif 'low' in file:
                    phenotype = 'Lower aggregation propensity in NKM'

                if 'BP' in file :
                    GO_type = 'Biological Process'
                elif 'CC' in file :
                    GO_type = 'Cellular Component'
                elif 'MF' in file :
                    GO_type = 'Molecular Function'
                if len(contents['overrepresentation']['group']) == 1:
                    try: 
                        GO = contents['overrepresentation']['group']['result']
                        if GO['input_list']['number_in_list'] < 2 :
                            print(f'{GO["result"]["term"]["label"]} excluded, {GO_type}')
                            pass
                        else: 
                            level = GO['term']['level']
                            go_id = GO['term']['id']
                            label = GO['term']['label']
                            fd = GO['input_list']['fold_enrichment']
                            pval = GO['input_list']['pValue']
                            protein_list = GO['input_list']['mapped_id_list']['mapped_id']
                            if protein_list is not str :
                                for ID in list(protein_list):
                                        all_gene_list.append([GO_type, level, go_id, label, ID, fd, pval, phenotype, analysis])
                            else: 
                                all_gene_list.append([GO_type, level, go_id, label, protein_list, fd, pval, phenotype, analysis])
                    except:
                        for lst in GO:
                            level = lst['term']['level']
                            go_id = lst['term']['id']
                            label = lst['term']['label']
                            fd = lst['input_list']['fold_enrichment']
                            pval = lst['input_list']['pValue']
                            protein_list = lst['input_list']['mapped_id_list']['mapped_id']
                            if protein_list is not str :
                                for ID in protein_list:
                                    all_gene_list.append([GO_type, level, go_id, label, ID, fd, pval, phenotype, analysis])
                            else: 
                                all_gene_list.append([GO_type, level, go_id, label, protein_list, fd, pval, phenotype, analysis])

                for GO in contents['overrepresentation']['group']:
                    if type(GO) is not str:
                        try:
                            if GO['result']['input_list']['number_in_list'] < 2 :
                                print(f'{GO["result"]["term"]["label"]} excluded, {GO_type}')
                                pass
                            else: 
                                level = GO['result']['term']['level']
                                go_id = GO['result']['term']['id']
                                label = GO['result']['term']['label']
                                fd = GO['result']['input_list']['fold_enrichment']
                                pval = GO['result']['input_list']['pValue']
                                protein_list = GO['result']['input_list']['mapped_id_list']['mapped_id']
                                if label in ['CCR chemokine receptor binding'] :
                                    for ID in list(protein_list):
                                        all_gene_list.append([GO_type, level, go_id, label, ID, fd, pval, phenotype, analysis])
                                if protein_list is not str :
                                    for ID in list(protein_list):
                                            all_gene_list.append([GO_type, level, go_id, label, ID, fd, pval, phenotype, analysis])
                                else: 
                                    all_gene_list.append([GO_type, level, go_id, label, protein_list, fd, pval, phenotype, analysis])
                        except:
                            try:
                                for lst in GO['result']:
                                    if lst['input_list']['number_in_list'] < 2 :
                                        print(f'{lst["term"]["label"]} excluded, {GO_type}')
                                        pass
                                    else: 
                                        level = lst['term']['level']
                                        go_id = lst['term']['id']
                                        label = lst['term']['label']
                                        fd = lst['input_list']['fold_enrichment']
                                        pval = lst['input_list']['pValue']
                                        protein_list = lst['input_list']['mapped_id_list']['mapped_id']
                                        if protein_list is not str :
                                            for ID in protein_list:
                                                all_gene_list.append([GO_type, level, go_id, label, ID, fd, pval, phenotype, analysis])
                                        else: 
                                            all_gene_list.append([GO_type, level, go_id, label, protein_list, fd, pval, phenotype, analysis])
                            except:
                                for res in GO['result']:
                                    if res['input_list']['number_in_list'] < 2 :
                                        print(f'{res["term"]["label"]} excluded, {GO_type}')
                                        pass

    all_GO_gene_list = pd.DataFrame(all_gene_list, columns=['GO Type', 'Level', 'GO ID', 'GO Term', 'proteinID', 'Fold Enrichment', 'raw p-value', 'Phenotype', 'Analysis'])
    all_GO_gene_list['description'] = all_GO_gene_list['proteinID'].progress_apply(check_annotations, args=(MM_fasta,))

    all_GO_gene_list['Subset'] = [ 'chaperone clients' if ID in mm_chap_clt.values else 'other proteins' for ID in all_GO_gene_list['proteinID'].values ]
    all_GO_gene_list['log2 Fold Enrichment'] = np.log2(all_GO_gene_list['Fold Enrichment'])
    all_GO_gene_list['-log10 p-value'] = -np.log10(all_GO_gene_list['raw p-value'])
    all_GO_gene_list.to_csv('../data/over-representation_analysis/stats/corrected_hypergeometric_tests/GO_gene_list.csv', index=False)
    return all_GO_gene_list


def multiple_chisquare_tests(terms_with_all_genes):
        all_tests = []
        for GO in np.unique(terms_with_all_genes['GO Term']):
            IN = terms_with_all_genes[terms_with_all_genes['GO Term'] == GO]
            CHAP_IN = IN[IN['Subset'] == 'chaperone clients']
            OTHERS_IN = IN[IN['Subset'] == 'other proteins']
            OUT = terms_with_all_genes[terms_with_all_genes['GO Term'] != GO]
            CHAP_OUT = OUT[OUT['Subset'] == 'chaperone clients']
            OTHERS_OUT = OUT[OUT['Subset'] == 'other proteins']
            OBSERVED = np.array([[len(CHAP_IN)/len(IN)*100, len(CHAP_OUT)/len(OUT)*100], [len(OTHERS_IN)/len(IN)*100, len(OTHERS_OUT)/len(OUT)*100]])
            EXPECTED = stats.contingency.expected_freq(OBSERVED)
            if (len(OBSERVED[OBSERVED >= 5]) > 0) & (len(EXPECTED[EXPECTED >= 5]) > 0):
                test_type = 'chisquare'
                chi2, pval, dof, expected = stats.chi2_contingency(OBSERVED)
                all_tests.append([GO, test_type, chi2, pval])
            else:
                test_type = 'barnard'
                barnard, pval = stats.barnard_exact(OBSERVED)
                all_tests.append([GO, test_type, barnard, pval])

        all_chisquare_pvals = pd.DataFrame(all_tests, columns=['GO Term', 'Test', 'statistic', 'p-value']).sort_values('p-value')
        all_chisquare_pvals['FDR'] = multipletests(all_chisquare_pvals['p-value'], alpha=0.05, method='fdr_bh')[1]
        all_chisquare_pvals.sort_values('FDR').to_csv('../data/over-representation_analysis/stats/corrected_chisquare/GO_distribution_chap_vs_others.csv', index=False)
        return all_chisquare_pvals

if __name__ == "__main__":
    tqdm.pandas()
    if not os.path.isfile('../data/over-representation_analysis/stats/corrected_hypergeometric_tests/GO_gene_list.csv'):
        mm_chap_clt = chaperone_clients_subset()
        terms_with_all_genes = build_gene_list(mm_chap_clt)
        print('Table S5 generated')
    else: 
        print('Table S5 already generated!')    
        terms_with_all_genes = pd.read_csv('../data/over-representation_analysis/stats/corrected_hypergeometric_tests/GO_gene_list.csv')

    GO = generate_input_for_figure()
    terms_with_all_genes = terms_with_all_genes[terms_with_all_genes['GO Term'].isin(GO['GO Term'])]
    hierarchical_order = list(pd.unique(terms_with_all_genes.sort_values(['GO Type', 'Fold Enrichment'], ascending=True)['GO Term']))[::-1]

    analysis_list = terms_with_all_genes['Analysis'].value_counts().index.tolist()
    analysis_cat = pd.Categorical(terms_with_all_genes['Analysis'], categories=analysis_list)
    terms_with_all_genes = terms_with_all_genes.assign(Analysis=analysis_cat)
    terms_with_all_genes['Analysis'] = terms_with_all_genes['Analysis'].cat.reorder_categories(['Domains', 'Proteins'])

    all_chisquare_pvals = multiple_chisquare_tests(terms_with_all_genes)
    print('GO Terms with significant chisquare pvalues ')
    print(all_chisquare_pvals[all_chisquare_pvals['FDR'] <= 0.05])

    fig = (p9.ggplot(
        terms_with_all_genes,
        p9.aes(x='GO Term', fill='Phenotype', alpha='Analysis')
        )
        + p9.geom_bar()
        + p9.scale_x_discrete(limits=hierarchical_order)
        + p9.scale_fill_manual(values=('red', 'blue'))
        + p9.labs(y='Protein count')
        + p9.guides(
            fill = p9.guide_legend(ncol=1),
            alpha = p9.guide_legend(ncol=1)
        )
        + p9.coord_flip()
        + p9.theme_classic()
        + p9.theme(figure_size=(6,15), 
                legend_background=p9.element_rect(size=2),
                legend_box='horizontal', 
                legend_position='top')
        + p9.facet_wrap('Subset')
    )

    fig.save('../figures/FIGURE_S1.png', dpi=300)
    # fig.save('../figures/FIGURE_S1.svg', dpi=300)
    # fig.save('../figures/FIGURE_S1.pdf', dpi=300)

    IN = terms_with_all_genes[terms_with_all_genes['GO Term'] == 'cytokine activity']
    CHAP_IN = IN[IN['Subset'] == 'chaperone clients']
    OTHERS_IN = IN[IN['Subset'] == 'other proteins']
    OUT = terms_with_all_genes[terms_with_all_genes['GO Term'] != 'cytokine activity']
    CHAP_OUT = OUT[OUT['Subset'] == 'chaperone clients']
    OTHERS_OUT = OUT[OUT['Subset'] == 'other proteins']
    OBSERVED = np.array([[len(CHAP_IN)/len(IN)*100, len(CHAP_OUT)/len(OUT)*100], [len(OTHERS_IN)/len(IN)*100, len(OTHERS_OUT)/len(OUT)*100]])
    EXPECTED = stats.contingency.expected_freq(OBSERVED)
    print('Observed')
    print(pd.DataFrame(OBSERVED, columns=['Chap', 'Others']))
    print('Expected')
    print(pd.DataFrame(EXPECTED, columns=['Chap', 'Others']))
    if (len(OBSERVED[OBSERVED >= 5]) > 0) & (len(EXPECTED[EXPECTED >= 5]) > 0):
        test_type = 'chisquare'
        chi2, pval, dof, expected = stats.chi2_contingency(OBSERVED)
        print([GO, test_type, chi2, pval])
    else:
        test_type = 'barnard'
        barnard, pval = stats.barnard_exact(OBSERVED)
        print.append([GO, test_type, barnard, pval])