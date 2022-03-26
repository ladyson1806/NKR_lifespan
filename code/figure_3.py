import os, re 
import numpy as np
import pandas as pd
import seaborn as sns
import plotnine as p9

def format_GO(x, format_type):
    if 'name' in format_type :
        GO_id = re.findall(r'( \(GO:[0-9]+\)$)', x)[0]
        x = x.replace(GO_id, '')
        return x
    if 'id' in format_type:
        GO_id = re.findall(r'\((GO:[0-9]+)\)$', x)[0]
        return GO_id


def build_master_table(analysis_type, agg_type):
    if 'domain' in analysis_type: 
        result_path = os.path.join('../data/over-representation_analysis/domains/OUTPUTS/')
    elif 'whole_sequence' in analysis_type:
        result_path = os.path.join(
            '../data/over-representation_analysis/whole_sequence/OUTPUTS/')

    if 'high' in agg_type :
        phenotype = 'Higher aggregation propensity in NKM'
    elif 'low' in agg_type:
        phenotype = 'Lower aggregation propensity in NKM'

    for file in os.listdir(result_path) :
        if agg_type in file :

            if 'CC' in file :
                CC_table = pd.read_table(os.path.join(result_path, file), comment='#', names=['Cellular Component', '# Reference list', '# Analyzed List', 'Expected', 'Over/Under', 'Fold Enrichment', 'raw p-value', 'FDR'])
                CC_table['GO_Type'] = ['Cellular Component' for i in CC_table.index]
                CC_table = CC_table.rename(columns={'Cellular Component': 'Gene Ontology Term'})
            elif 'MF' in file :
                MF_table = pd.read_table(os.path.join(result_path, file), comment='#', names=['Molecular Function', '# Reference list', '# Analyzed List', 'Expected', 'Over/Under', 'Fold Enrichment', 'raw p-value', 'FDR'])
                MF_table['GO_Type'] = ['Molecular Function' for i in MF_table.index]
                MF_table = MF_table.rename(columns={'Molecular Function': 'Gene Ontology Term'})
            elif 'BP' in file :
                BP_table = pd.read_table(os.path.join(result_path, file), comment='#', names=['Biological Process', '# Reference list', '# Analyzed List', 'Expected', 'Over/Under', 'Fold Enrichment', 'raw p-value', 'FDR'])
                BP_table['GO_Type'] = ['Biological Process' for i in BP_table.index]
                BP_table = BP_table.rename(columns={'Biological Process': 'Gene Ontology Term'})

    all_GO_table = pd.concat([CC_table, MF_table, BP_table])
    all_GO_hits = all_GO_table[(all_GO_table['FDR'] <= 0.05) & (all_GO_table['Fold Enrichment'] != ' < 0.01')].reset_index(drop=True)
    all_GO_hits['Fold Enrichment'] = all_GO_hits['Fold Enrichment'].astype(float)
    all_GO_hits['Phenotype'] = phenotype 
    all_GO_hits = all_GO_hits.sort_values(['Fold Enrichment'], ascending=False)
    return all_GO_hits


def nr_GO_terms():
    all_revigo = []
    revigo_path = '../data/over-representation_analysis/revigo'
    revigo_columns = ['TermID',  'Name',  'Frequency',  'PlotX',  'PlotY',  'LogSize',  'Value',  'Uniqueness',  'Dispensability',  'Representative',  'Eliminated']
    for file in os.listdir(revigo_path):
        if 'inputs' not in file:
            all_revigo.append(pd.read_csv(os.path.join(revigo_path, file), names=revigo_columns)[1:])

    REVIGO = pd.concat(all_revigo).reset_index(drop=True)
    return list(REVIGO[REVIGO['Eliminated'].str.contains('False')]['TermID'])


def generate_input_for_figure():
    HD = build_master_table('domain', 'high')
    LD = build_master_table('domain', 'low')
    HWS = build_master_table('whole_sequence', 'high')
    LWS = build_master_table('whole_sequence', 'low')

    #### Domains
    ALL_DOM = pd.concat([HD, LD]).reset_index(drop=True).sort_values('GO_Type')
    ALL_DOM['log2 Fold Enrichment'] = np.log2(ALL_DOM['Fold Enrichment'])
    ALL_DOM['-log10(p-value)'] = -np.log10(ALL_DOM['raw p-value'])
    ALL_DOM['GO Term'] = ALL_DOM['Gene Ontology Term'].apply(format_GO, args=('name',))
    ALL_DOM['GO ID'] = ALL_DOM['Gene Ontology Term'].apply(format_GO, args=('id',))
    ALL_DOM['Analysis'] = 'Domains'
    ALL_DOM = ALL_DOM.sort_values(['Phenotype', 'GO_Type', 'Fold Enrichment'], ascending=False).reset_index(drop=True)
    ALL_DOM[['GO Term', 'GO ID', 'GO_Type', 'Phenotype', '# Reference list', '# Analyzed List', 'Expected', 'Over/Under', 'Fold Enrichment', 'raw p-value', 'FDR']].to_csv('../data/over-representation_analysis/stats/corrected_hypergeometric_tests/all_GO_domain_agg.csv', index=False)
    print('Table S2 generated')
    #### Proteins
    ALL_PROT = pd.concat([HWS, LWS]).reset_index(drop=True)
    ALL_PROT['log2 Fold Enrichment'] = np.log2(ALL_PROT['Fold Enrichment'])
    ALL_PROT['-log10(p-value)'] = -np.log10(ALL_PROT['raw p-value'])
    ALL_PROT['GO Term'] = ALL_PROT['Gene Ontology Term'].apply(format_GO, args=('name',))
    ALL_PROT['GO ID'] = ALL_PROT['Gene Ontology Term'].apply(format_GO, args=('id',))
    ALL_PROT['Analysis'] = 'Proteins'
    ALL_PROT = ALL_PROT.sort_values(['Phenotype', 'GO_Type', 'Fold Enrichment'], ascending=False).reset_index(drop=True)
    ALL_PROT[['GO Term', 'GO ID', 'GO_Type', 'Phenotype', '# Reference list', '# Analyzed List', 'Expected', 'Over/Under', 'Fold Enrichment', 'raw p-value', 'FDR']].to_csv('../data/over-representation_analysis/stats/corrected_hypergeometric_tests/all_GO_protein_agg.csv', index=False)
    print('Table S3 generated')

    ALL = pd.concat([ALL_DOM, ALL_PROT])

    #### Removing terms with less than 5 proteins
    ALL = ALL[ALL['# Analyzed List'] >= 5]
    #### Removing redundant terms with Revigo
    # ALL = ALL[ALL['GO ID'].isin(revigo_terms)] 

    ### Panel order
    analysis_list = ALL['Analysis'].value_counts().index.tolist()
    analysis_cat = pd.Categorical(ALL['Analysis'], categories=analysis_list)
    ALL = ALL.assign(Analysis = analysis_cat)
    ALL['Analysis'] = ALL['Analysis'].cat.reorder_categories(['Proteins', 'Domains'])

    return ALL.sort_values(['GO_Type', 'Fold Enrichment','Analysis', 'Phenotype'])


if __name__ == "__main__":
    terms_with_genes = pd.read_csv('../data/over-representation_analysis/GO_gene_list.csv')

    GO = generate_input_for_figure()
    hierarchical_order = list(pd.unique(GO[GO['GO ID'].isin(terms_with_genes['GO ID'])]['GO Term']))
    GO_with_genes = GO[GO['GO Term'].isin(hierarchical_order)].sort_values(['GO_Type', 'Fold Enrichment'], ascending=True)[::-1]

    fig = (p9.ggplot(GO_with_genes, p9.aes(x='log2 Fold Enrichment', y='GO Term', size='-log10(p-value)', color='Phenotype', shape='Analysis'))
        + p9.geom_point()
        + p9.scale_x_continuous(breaks=np.arange(-4, 4.5, 1), limits=[-4, 4.5])
        + p9.scale_y_discrete(limits=pd.unique(GO_with_genes['GO Term']))
        + p9.scale_color_manual(values=("red", "blue"))
        + p9.guides(
            shape=p9.guide_legend(ncol=1, title='Analysis type'),
            color=p9.guide_legend(ncol=1, title='Phenotype'),
            size=p9.guide_legend(ncol=1, title='-log10(p-value)')
    )
        + p9.geom_vline(
            xintercept=0,
            linetype="dashed"
    )
        + p9.labs(
            x='',
            y="GO Terms",
    )
        + p9.theme_classic()
        + p9.theme(
            figure_size=(6,15), 
            legend_position='top',
            legend_box='horizontal',
            panel_grid_major_x=p9.element_line(color='white')
        )
    )
    
    fig.save('../figures/FIGURE_3.png', dpi=300)
    # fig.save('../figures/FIGURE_3.svg', dpi=300)
    # fig.save('../figures/FIGURE_3.pdf', dpi=300)