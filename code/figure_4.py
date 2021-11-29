import os, progressbar, re, subprocess

import numpy as np
import pandas as pd
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import statsmodels.sandbox.stats.multicomp as multicomp


def is_significant(mut_tol):
    if mut_tol > 2 :
        return 'HIGH'
    if mut_tol < -2:
        return 'LOW'
    else :
        return 'non significant'


def cluster(x):
    if x in HG_high_agg['proteinID_x'].values:
        return 'HIGH'
    else :
        return 'LOW'


def generate_figure_4A(prot_muttol_table):
    sns.set_context("paper", font_scale=2)
    sns.set_style("ticks") 
    sns.despine(offset=20)

    fig = plt.figure(figsize= (8, 8))
    sns.scatterplot(x=prot_muttol_table['mutTol_x'], y=prot_muttol_table['mutTol_y'], hue=prot_muttol_table['MT_DIFF'], hue_order=['non significant', 'HIGH', 'LOW'], palette=['black', 'purple', 'salmon'],legend=False)
    plt.ylabel('Mouse mutation tolerance')
    plt.xlabel('Naked mole-rat mutation tolerance')
    
    fig.savefig('../figures/FIGURE4_A.png', format='png', dpi=300)
    # fig.savefig('../figures/FIGURE4_A.svg', format='svg', dpi=300)
    # fig.savefig('../figures/FIGURE4_A.pdf', format='pdf', dpi=300)

    print('Stats for Figure 4A')
    corr, pval = stats.pearsonr(x=prot_muttol_table['mutTol_x'], y=prot_muttol_table['mutTol_y'])
    print('Correlation between mutation tolerance')
    print(f'correlation score:{corr}, pvalue:{pval} \n')


def generate_figure_4BC(prot_muttol_table):
    mm_mutTol = prot_muttol_table[ (prot_muttol_table['diff_mut_z-scores'] < -2)]
    hg_mutTol = prot_muttol_table[ (prot_muttol_table['diff_mut_z-scores'] > 2)]

    plt.rcParams["figure.figsize"] = (18,9)
    fig, axes = plt.subplots(nrows=1, ncols=2, sharex=False, sharey=True)
    sns.set_context("paper", font_scale=2)
    sns.set_style("ticks") 
    sns.despine(offset=10)

    sns.regplot(x=prot_muttol_table['Aggregation_y'], y=prot_muttol_table['mutTol_y'], color='black', ax=axes[0])
    sns.regplot(x=mm_mutTol['Aggregation_y'], y=mm_mutTol['mutTol_y'], color='salmon', ax=axes[0])
    axes[0].set_xlabel('Mouse whole-protein sequence \n aggregation propensity score')
    axes[0].set_ylabel('Mouse mutation tolerance')
    axes[0].set_xlim([0, 25])
    axes[0].set_ylim([0, 1])

    sns.regplot(x=prot_muttol_table['Aggregation_x'], y=prot_muttol_table['mutTol_x'], color='black', ax=axes[1])
    sns.regplot(x=hg_mutTol['Aggregation_x'], y=hg_mutTol['mutTol_x'], color='purple')
    plt.xlabel('Naked mole-rat whole-protein sequence \n aggregation propensity score', )
    plt.ylabel('Naked mole-rat mutation tolerance', labelpad=5)
    axes[1].set_xlim([0, 25])
    axes[1].set_ylim([0, 1])
    axes[1].tick_params(axis='both', direction='out', length=6, width=1, colors='black', grid_color='black', grid_alpha=0.5, labelleft=True)

    fig.savefig('../figures/FIGURE4_BC.png', format='png', dpi=300)
    # fig.savefig('../figures/FIGURE4_BC.svg', format='svg', dpi=300)
    # fig.savefig('../figures/FIGURE4_BC.pdf', format='pdf', dpi=300)

    print('\nStats for Figure 4BC')
    print('Correlation between whole-protein sequence aggregation propensity and mutation tolerance in mouse')
    corr, pval = stats.pearsonr(x=prot_muttol_table['mutTol_y'], y=prot_muttol_table['Aggregation_y'])
    print(f'correlation score:{corr}, pvalue:{pval} (M) \n')

    print('Difference of distribution in aggregation propensity between proteins with significant difference of mutation tolerance and the rest of the proteins in mouse')
    print('KS test:', stats.kstest(mm_mutTol['Aggregation_y'], prot_muttol_table['Aggregation_y']), '\n')

    print('Correlation between whole-protein sequence aggregation propensity and mutation tolerance in naked mole-rat')
    corr, pval = stats.pearsonr(x=prot_muttol_table['mutTol_x'], y=prot_muttol_table['Aggregation_x'])
    print(f'correlation score:{corr}, pvalue:{pval} (NKM)\n')

    print('Difference of distribution in  aggregation propensity between proteins with significant difference of mutation tolerance and the rest of the proteins in naked mole-rat')
    print('KS test:', stats.kstest(hg_mutTol['Aggregation_x'], prot_muttol_table['Aggregation_x']), '\n')


def generate_figure_4D(sign_agg_table, prot_muttol_table):
    sns.set_context("paper", font_scale=2)
    sns.set_style("ticks") 
    sns.despine(offset=20)

    plt.rcParams["figure.figsize"] = (12,12)

    j = sns.jointplot(data=sign_agg_table, y='mutTol_y', x='Aggregation_y', hue='AGG_DIFF', palette=['red', 'blue'], legend=False)
    j.set_axis_labels('Mouse whole-protein sequence \naggregation propensity score', 'Mouse mutation tolerance', fontsize=16)
    j.savefig('../figures/FIGURE4_D.png', format='png', dpi=300)
    # j.savefig('../figures/FIGURED_E.svg', format='svg', dpi=300)
    # j.savefig('../figures/FIGURED_E.pdf', format='pdf', dpi=300)

    MT_vs_HIGH_AGG = prot_muttol_table[prot_muttol_table['proteinID_x'].isin(HG_high_agg['proteinID_x'])]
    MT_vs_LOW_AGG = prot_muttol_table[prot_muttol_table['proteinID_x'].isin(HG_low_agg['proteinID_x'])]

    print('\nStats for Figure 4D')
    print('Difference of distribution in mutation tolerance between proteins with sgnificant difference of aggregation propensity in mouse')
    print('KS test:', stats.kstest(MT_vs_HIGH_AGG['mutTol_y'], MT_vs_LOW_AGG['mutTol_y']), '\n')


def generate_figure_4E(sign_agg_table, prot_muttol_table):
    sns.set_context("paper", font_scale=2)
    sns.set_style("ticks") 
    sns.despine(offset=20)

    plt.rcParams["figure.figsize"] = (12,12)

    j = sns.jointplot(data=sign_agg_table, y='mutTol_x', x='Aggregation_x', hue='AGG_DIFF', palette=['red', 'blue'], legend=False)
    j.set_axis_labels('Naked mole-rat whole-protein sequence \naggregation propensity score', 'Naked mole-rat mutation tolerance', fontsize=16)
    j.savefig('../figures/FIGURE4_E.png', format='png', dpi=300)
    # j.savefig('../figures/FIGURE4_E.svg', format='svg', dpi=300)
    # j.savefig('../figures/FIGURE4_E.pdf', format='pdf', dpi=300)

    MT_vs_HIGH_AGG = prot_muttol_table[prot_muttol_table['proteinID_x'].isin(HG_high_agg['proteinID_x'])]
    MT_vs_LOW_AGG = prot_muttol_table[prot_muttol_table['proteinID_x'].isin(HG_low_agg['proteinID_x'])]

    print('\nStats for Figure 4E')
    print('Difference of distribution in mutation tolerance between proteins with sgnificant difference of aggregation propensity in naked mole-rat')
    print('KS test:',stats.kstest(MT_vs_HIGH_AGG['mutTol_x'], MT_vs_LOW_AGG['mutTol_x'], '\n'))

if __name__ == "__main__":
    #### Table with all per-protein aggregation propensity scores
    prot_agg_table = pd.read_csv('../data/aggregation_propensity/HGMM_agg_scores.csv', sep=',')
    prot_agg_table['delta_aggregation'] = prot_agg_table['Aggregation_x'] - prot_agg_table['Aggregation_y']
    prot_agg_table['delta_agg_z-scores'] = stats.zscore(prot_agg_table['delta_aggregation'])
    prot_agg_table['difference_of_aggregation'] = prot_agg_table['delta_agg_z-scores'].apply(is_significant)

    prot_muttol_table = pd.read_csv('../data/mutation_tolerance/all_mt_scores.csv', sep='\t')
    prot_muttol_table['diff_mut'] = prot_muttol_table['mutTol_x'] - prot_muttol_table['mutTol_y']
    prot_muttol_table['diff_mut_z-scores'] = stats.zscore(prot_muttol_table['diff_mut'])
    prot_muttol_table['MT_DIFF'] = prot_muttol_table['diff_mut_z-scores'].apply(is_significant)
    
    HG_high_agg = prot_agg_table[prot_agg_table['delta_agg_z-scores'] > 2]
    HG_low_agg = prot_agg_table[prot_agg_table['delta_agg_z-scores'] < -2]

    sign_agg_table = prot_muttol_table[(prot_muttol_table['proteinID_x'].isin(HG_high_agg['proteinID_x'])) | (prot_muttol_table['proteinID_x'].isin(HG_low_agg['proteinID_x'])) ].reset_index(drop=True)
    sign_agg_table['AGG_DIFF'] = sign_agg_table['proteinID_x'].apply(cluster)

    generate_figure_4A(prot_muttol_table)
    generate_figure_4BC(prot_muttol_table)
    generate_figure_4D(sign_agg_table, prot_muttol_table)
    generate_figure_4E(sign_agg_table, prot_muttol_table)