import pandas as pd
import scanpy as sc
from anndata import AnnData


def _check_expression_mtx(adata):
    m = adata.X < 0
    if m.size > 0:
        return False


def get_annotated(adata, genes):
    if _check_expression_mtx(adata):
        print('Expression matrix seems to be Z transformed!!!')
    else:
        if isinstance(genes, list):
            # Add expression levels to df
            subdata = adata[:, adata.var_names.isin(genes)].copy()
            df = subdata.to_df()
            df = pd.concat([subdata.obs, df], axis=1)

            # Annotate Positivity
            for gene in genes:
                if gene == 'MS4A1':
                    df['CD20_status'] = df.MS4A1.apply(lambda x: 'positive' if x > 0 else 'negative')

                else:
                    gname = f'{gene}_status'
                    df[gname] = df[gene].apply(lambda x: 'positive' if x > 0 else 'negative')

            # Add annotation to adata
            for gene in genes:
                if gene == 'MS4A1':
                    adata.obs['CD20_status'] = df['CD20_status'].copy()

                else:
                    gname = f'{gene}_status'
                    adata.obs[gname] = df[gname].copy()

            return df, adata


def pct_table(df, group, file_name='table.xlsx', gene='MS4A1', threshold=0):
    from collections import OrderedDict
    vals = OrderedDict({})
    for cell in df[group].value_counts().sort_index().index:
        allcells = df[df[group] == cell].shape[0]
        cds = df[(df[group] == cell) & (df[gene] > threshold)].shape[0]
        if allcells == 0:
            pct = 0
        else:
            pct = round(cds / allcells * 100, 2)
        vals[cell] = [allcells, cds, pct]
    table = pd.DataFrame(vals, index=['allcells', f'{gene} Expressing', 'Percentage'])
    table.to_excel(file_name)
    return table


def contig(df, group1, group2, file_name='table.xlsx', gene='MS4A1', threshold=0):
    from collections import OrderedDict
    grp = OrderedDict({})
    for i in df[group1].value_counts().sort_index().index:
        grp[i] = pct_table(df[df[group1] == i], group2, file_name=file_name, gene=gene, threshold=threshold)
    table = pd.concat(grp)
    table.to_excel(file_name)
    return table


def analyze_pct(df, groups, genes, threshold=0, folder_name='pct'):
    from itertools import combinations
    import os
    import shutil

    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)

    # Calculate percentage positives for groupby in obs
    for gen in genes:
        for grp in groups:
            pct_table(df=df,
                      group=grp,
                      file_name=f'{folder_name}/{grp}-{gen}.xlsx',
                      gene=gen,
                      threshold=threshold
                      )
    for gen in genes:
        for g1, g2 in combinations(groups, 2):
            contig(df=df,
                   group1=g1,
                   group2=g2,
                   file_name=f'{folder_name}/{g1}+{g2}-{gen}.xlsx',
                   gene=gen,
                   threshold=threshold
                   )

    shutil.make_archive(folder_name, 'zip', folder_name)


def analyze_dge(adata: AnnData, label_keys: list, factors: list, versus: list, analyze_global: bool = True,
                folder_name: str = 'dge'):
    import os
    import shutil

    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)

# Iterate over labels
    bdata = adata.copy()
    for label_key in label_keys:
        for label in adata.obs[label_key].unique():
            adata = bdata.copy()
            adata = adata[adata.obs[label_key] == label].copy()
            if analyze_global:
                for vs in versus:
                    sc.tl.rank_genes_groups(adata, vs, method='wilcoxon')
                    dge = sc.get.rank_genes_groups_df(adata, 'positive', pval_cutoff=0.05).sort_values('logfoldchanges',
                                                                                                       ascending=False)
                    dge.to_excel(f"global-{label_key}={label}-{vs.split('_', 1)[0]}+_vs_{vs.split('_', 1)[0]}-.xlsx")

            cdata = adata.copy()
            for factor in factors:

                for lev in adata.obs[factor].unique():
                    adata = cdata.copy()
                    adata = adata[adata.obs[factor] == lev].copy()

                    for vs in versus:
                        sc.tl.rank_genes_groups(adata, vs, method='wilcoxon')
                        dge = sc.get.rank_genes_groups_df(adata, 'positive',
                                                          pval_cutoff=0.05).sort_values('logfoldchanges',
                                                                                        ascending=False)
                        dge.to_excel(
                            f"{folder_name}/{label_key}={label}+{factor}={lev}-{vs.split('_', 1)[0]}+_vs_{vs.split('_', 1)[0]}-.xlsx")

    shutil.make_archive(folder_name, 'zip', folder_name)
