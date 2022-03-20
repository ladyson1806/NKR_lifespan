import numpy as np
import pandas as pd 
import seaborn as sns
import scipy.stats as stats
import matplotlib.pyplot as plt
import ete3, os, subprocess

from ete3 import Tree, faces, TreeStyle


def group_color(name):
    name = name.lower() 
    if name in ["heterocephalus_glaber", "fukomys_damarensis", "cavia_porcellus", 'octodon_degus', 'chinchilla_lanigera']:
        return "blue"
    elif name in ["rattus_norvegicus", "mus_musculus", "meriones_unguiculatus", "mesocricetus_auratus", "microtus_ochrogaster", "peromyscus_maniculatus", "neotoma_lepida", "jaculus_jaculus"]:
        return "red"
    elif name in ["dipodomys_ordii", "castor_canadensis"]:
        return "green"
    elif name in ["ictidomys_tridecemlineatus", "marmota_flaviventris", "marmota_marmota"]:
        return "orange"

def generate_tree(rodent_table, rodent_tree):
    def get_lgv(name):
        return rodent_table[rodent_table['Scientific name'] == name]['Maximum longevity (yrs)'].values[0]

    def my_layout(node):
        if node.is_leaf():
            ## Get lifespan + family color information
            longevity = get_lgv(node.name.lower())
            fam_color = group_color(node.name)
            
            ## Name formating
            if node.name.lower() in ['mus_musculus', 'heterocephalus_glaber']:
                longNameFace = faces.TextFace(node.name.replace('_', ' '), ftype='Arial', fsize=10, fgcolor=fam_color, fstyle='italic', bold=True)
                longNameFace.margin_left = 10
            else: 
                longNameFace = faces.TextFace(node.name.replace('_', ' '), ftype='Arial', fsize=10, fgcolor=fam_color, fstyle='italic')
                longNameFace.margin_left = 10
            faces.add_face_to_node(longNameFace, node, column=0)
            
            ## Longevity display
            if node.name.lower() in ['mus_musculus', 'heterocephalus_glaber']:
                L = faces.TextFace(f'{longevity} years', ftype='Arial', fsize=10, bold=True)
            else:
                L = faces.TextFace(f'{longevity} years', ftype='Arial', fsize=10)
            L.margin_left = 20
            L.margin_top= 2
            faces.add_face_to_node(L, node, column=1)
    
        nstyle = ete3.NodeStyle()
        nstyle["size"] = 0
        nstyle["fgcolor"] = 'white'
        node.set_style(nstyle)

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.layout_fn = my_layout

    rodent_tree.render("%%inline", w=800, tree_style=ts)
    rodent_tree.render("../figures/FIGURE1_A_v2.png", dpi=300, tree_style=ts)
    # rodent_tree.render("../figures/FIGURE1_A.svg", dpi=300, tree_style=ts)
    # rodent_tree.render("../figures/FIGURE1_A.pdf", dpi=300, tree_style=ts)

def generate_scatterplots(rodent_table):
    sns.set_context("paper", font_scale=2)
    sns.set_style("ticks") 
    sns.despine(offset=20)
    fig,axes = plt.subplots(nrows=1, ncols=3, figsize=(15, 5), sharey=True)

    focus = (rodent_table['Common name'] == 'House mouse') | (rodent_table['Common name'] == 'Naked mole-rat')
    rodent_table['Rodents'] = np.where(focus==True, "Chosen species", "Other rodents")
    rodent_table['color'] = rodent_table['Scientific name'].apply(group_color)
    markers = {"Chosen species": "s", "Other rodents": "X"}

    sns.scatterplot(y=np.log10(rodent_table['Maximum longevity (yrs)']), x=np.log10(rodent_table['Body mass (g)']), hue=rodent_table['color'], hue_order=['blue', 'red', 'green', 'orange'], palette=['blue', 'red', 'green', 'orange'], style=rodent_table['Rodents'], markers=['s','o'], style_order=['Chosen species', 'Other rodents'], legend=False, ax=axes[0])
    axes[0].set_xlabel('log10 Adult weight (gram)', fontsize=12)
    axes[0].set_ylabel('log10 Maximum lifespan (year)', fontsize=12, labelpad=15)


    sns.scatterplot(y=np.log10(rodent_table['Maximum longevity (yrs)']), x=np.log10(rodent_table['Metabolic rate (W)']), hue=rodent_table['color'], hue_order=['blue', 'red', 'green', 'orange'], palette=['blue', 'red', 'green', 'orange'], style=rodent_table['Rodents'], markers=['s','o'], style_order=['Chosen species', 'Other rodents'], legend=False, ax=axes[1])
    axes[1].set_xlabel('log10 Basal metabolic rate (W)', fontsize=12)
    axes[1].set_ylabel('log10 Maximum lifespan (year)', fontsize=12, labelpad=15)
    axes[1].tick_params(axis='both', direction='out', length=6, width=1, colors='black', grid_color='black', grid_alpha=0.5, labelleft=True)

    sns.scatterplot(y=np.log10(rodent_table['Maximum longevity (yrs)']), x=np.log10(rodent_table['Female maturity (days)']), hue=rodent_table['color'], hue_order=['blue', 'red', 'green', 'orange'], palette=['blue', 'red', 'green', 'orange'], style=rodent_table['Rodents'], markers=['s','o'], style_order=['Chosen species', 'Other rodents'], legend=False,ax=axes[2])
    axes[2].set_xlabel('log10 Female maturity (day)', fontsize=12)
    axes[2].set_ylabel('log10 Maximum lifespan (year)', fontsize=12, labelpad=15)
    axes[2].tick_params(axis='both', direction='out', length=6, width=1, colors='black', grid_color='black', grid_alpha=0.5, labelleft=True)

    fig.savefig('../figures/FIGURE1_BCD_v2.png', format='png')
    # fig.savefig('../figures/FIGURE1_BCD.svg', format='svg')
    # fig.savefig('../figures/FIGURE1_BCD.pdf', format='pdf')

    rodent_table = rodent_table.dropna()
    corr, pval = stats.pearsonr(x=rodent_table['Body mass (g)'], y=rodent_table['Maximum longevity (yrs)'])
    print(f'correlation score:{corr}, pvalue:{pval} (Body mass vs. Maximum lifespan )')

    corr, pval = stats.pearsonr(x=rodent_table['Female maturity (days)'], y=rodent_table['Maximum longevity (yrs)'])
    print(f'correlation score:{corr}, pvalue:{pval} (Female maturity vs. Maximum lifespan )')

    corr, pval = stats.pearsonr(x=rodent_table['Metabolic rate (W)'], y=rodent_table['Maximum longevity (yrs)'])
    print(f'correlation score:{corr}, pvalue:{pval} (Metabolic rate vs. Maximum lifespan )')


if __name__ == "__main__":
    rodent_tree = ete3.Tree('../data/phylogeny/rodent_in_anage_phylogeny.nwk', format=1)
    rodent_anage_table = pd.read_csv('../data/anage/rodent_table_v2.csv')
    #### Generate Figure 1
    generate_tree(rodent_anage_table, rodent_tree)
    #### Generate Figure 1 B,C,D
    generate_scatterplots(rodent_anage_table)