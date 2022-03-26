import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
import json, os, progressbar, re, time

from Bio import SeqIO
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from matplotlib_venn import venn3, venn3_circles
from matplotlib_venn import venn2, venn2_circles
from plotnine import *



def is_significant(agg_score):
    if agg_score > 2:
        return 'red'
    if agg_score < -2:
        return 'blue'
    else:
        return 'black'
  
def generate_figure(prot_agg_table, dom_agg_table, mm_chap_clt):
    sns.set_context("paper", font_scale=2)
    sns.set_style("ticks") 
    sns.despine(offset=20)

    fig,axes = plt.subplots(nrows=2, ncols=2, figsize=(20, 20), sharex=False, sharey=False)

    sns.scatterplot(x=prot_agg_table['Aggregation_x'], y=prot_agg_table['Aggregation_y'], hue=prot_agg_table['difference_of_aggregation'], palette=['black', 'blue', 'red'], alpha=0.75, ax=axes[0, 1], legend=False)
    axes[0,0].set_ylabel('Mouse whole-protein sequence \naggregation propensity score', labelpad=8)
    axes[0,0].set_xlabel('Naked-mole rat whole-protein sequence \naggregation propensity score', visible=True)
    axes[0,0].set_title('All proteins',  fontstyle='italic', loc='left')
    axes[0,0].set_xlim(0, 50)   
    axes[0,0].set_ylim(0, 50)  

    CHAP = prot_agg_table[prot_agg_table['proteinID_y'].isin(mm_chap_clt)].sort_values('difference_of_aggregation')
    sns.scatterplot(x=CHAP['Aggregation_x'], y=CHAP['Aggregation_y'], hue=CHAP['difference_of_aggregation'], palette=['black', 'blue', 'red'], alpha=0.75, ax=axes[1, 1], legend=False)
    axes[0,1].set_ylabel('Mouse whole-protein sequence \naggregation propensity score')
    axes[0,1].set_xlabel('Naked-mole rat whole-protein sequence \naggregation propensity score')
    axes[0,1].set_title('Chaperone client proteins', fontstyle='italic', loc='left')
    axes[0,1].set_xlim(0, 30)   
    axes[0,1].set_ylim(0, 30) 

    sns.scatterplot(y=dom_agg_table['dom_agg_score_y'], x=dom_agg_table['dom_agg_score_x'], hue=dom_agg_table['difference_of_aggregation'], palette=['black', 'blue', 'red'], ax=axes[0, 0], legend=False)
    axes[1,0].set_ylabel('Mouse per-domain \naggregation propensity score', labelpad=8)
    axes[1,0].set_xlabel('Naked-mole rat per-domain \naggregation propensity score', visible=True)
    axes[1,0].set_title('Domains in all proteins', fontstyle='italic', loc='left')
    axes[1,0].set_xlim(0, 40)   
    axes[1,0].set_ylim(0, 40) 

    CHAP_DOM = dom_agg_table[dom_agg_table['proteinID_y'].isin(mm_chap_clt)].sort_values('difference_of_aggregation')
    sns.scatterplot(y=CHAP_DOM['dom_agg_score_y'], x=CHAP_DOM['dom_agg_score_x'], hue=CHAP_DOM['difference_of_aggregation'], palette=['black', 'blue', 'red'], ax=axes[1, 0], legend=False)
    axes[1,1].set_ylabel('Mouse per-domain \naggregation propensity score', labelpad=8)
    axes[1,1].set_xlabel('Naked-mole rat per-domain \naggregation propensity score')
    axes[1,1].set_title('Domains in chaperone client proteins', fontstyle='italic', loc='left')
    axes[1,1].set_xlim(0, 15)   
    axes[1,1].set_ylim(0, 15) 
    
    fig.savefig('../figures/FIGURE_2.png', format='png', quality=300)
    # fig.savefig('../figures/FIGURE_2.svg', format='svg', quality=300)
    # fig.savefig('../figures/FIGURE_2.svg', format='svg', quality=300)


    #### Correlation 

    print('Correlation between whole-protein sequence aggregation propensity')
    ## Correlation between HG and MM Tango scores - All proteins
    corr, pval = stats.pearsonr(prot_agg_table['Aggregation_x'], prot_agg_table['Aggregation_y'])
    print(f'correlation score:{corr}, pvalue:{pval} (All dataset)')

    ## Correlation between HG and MM Tango scores - Chaperone client proteins
    corr, pval = stats.pearsonr(CHAP['Aggregation_x'], CHAP['Aggregation_y'])
    print(f'correlation score:{corr}, pvalue:{pval} (Chaperone client proteins)')


    print('\nCorrelation between per-domain aggregation propensity')
    ## Correlation between HG and MM Tango scores - All domains
    corr, pval = stats.pearsonr(dom_agg_table['dom_agg_score_x'], dom_agg_table['dom_agg_score_y'])
    print(f'correlation score:{corr}, pvalue:{pval} (All dataset)')

    ## Correlation between HG and MM Tango scores - Domains in chaperone client proteins
    corr, pval = stats.pearsonr(CHAP_DOM['dom_agg_score_x'], CHAP_DOM['dom_agg_score_y'])
    print(f'correlation score:{corr}, pvalue:{pval} (Chaperone client proteins)')

    print('\n')
    #### T-tests 

    print('Difference of delta agg distribution for whole-protein sequence scores in chaperone clients and the rest of the proteins')
    full_stat, full_pval = stats.ttest_ind(prot_agg_table[~prot_agg_table['proteinID_y'].isin(CHAP['proteinID_y'])]['delta_agg_z-scores'], CHAP['delta_agg_z-scores'])
    print(full_stat, full_pval)

    print('Difference of delta agg distribution for domain scores for chaperone clients and the rest of the proteins')
    dom_stat, dom_pval = stats.ttest_ind(dom_agg_table[~dom_agg_table['proteinID_y'].isin(CHAP_DOM['proteinID_y'])]['delta_dom_agg_z-scores'], CHAP_DOM['delta_dom_agg_z-scores'])
    print(dom_stat, dom_pval)

if __name__ == "__main__":
    #### Table with all per-protein aggregation propensity scores
    prot_agg_table = pd.read_csv('../data/aggregation_propensity/HGMM_agg_scores.csv', sep=',')
    prot_agg_table['delta_aggregation'] = prot_agg_table['Aggregation_x'] - prot_agg_table['Aggregation_y']
    prot_agg_table['delta_agg_z-scores'] = stats.zscore(prot_agg_table['delta_aggregation'])
    prot_agg_table['difference_of_aggregation'] = prot_agg_table['delta_agg_z-scores'].apply(is_significant)

    #### Table with all per-domain aggregation propensity scores
    dom_agg_table = pd.read_csv('../data/aggregation_propensity/HGMM_dom_agg_scores.csv', sep='\t')
    dom_agg_table['delta_dom_aggregation'] = dom_agg_table['dom_agg_score_x'] - dom_agg_table['dom_agg_score_y']
    dom_agg_table['delta_dom_agg_z-scores'] = stats.zscore(dom_agg_table['delta_dom_aggregation'])
    dom_agg_table['difference_of_aggregation'] = dom_agg_table['delta_dom_agg_z-scores'].apply(is_significant)


    #### List of chaperone client proteins
    uniprot_mapping = pd.read_csv('../data/chaperone_clients/human_ensembl_to_uniprot.tab', sep='\t')
    hs_mm_orthologs = pd.read_csv('../data/chaperone_clients/HS_MM_uni_ortholog_groups.csv', sep='\t')
    hs_mm_orthologs = hs_mm_orthologs[['proteinID_x', 'proteinID_y']]
    mm_chap_clt = hs_mm_orthologs[hs_mm_orthologs['proteinID_x'].isin(uniprot_mapping['Entry'])]['proteinID_y']

    generate_figure(prot_agg_table.sort_values('difference_of_aggregation'), dom_agg_table.sort_values('difference_of_aggregation'), mm_chap_clt)

