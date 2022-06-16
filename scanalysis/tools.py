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


def analyze_pct(df, label_keys: list, group_keys, genes: list, analyze_global: bool=True, threshold=0, folder_name='pct'):
    from itertools import combinations
    import os
    import shutil

    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)

    bdf = df.copy()
    for label_key in label_keys:
        for gene in genes:
            if analyze_global:
                pct_table(df,
                          group=label_key,
                          gene=gene,
                          file_name=f"{folder_name}/global-{label_key}-{gene}.xlsx")

        for label in df[label_key].unique():
            df = bdf.copy()
            df = df[df[label_key] == label].copy()

            for gene in genes:

                for group_key in group_keys:
                    pct_table(df,
                              group=group_key,
                              file_name=f"{folder_name}/{label_key}={label}--{group_key}-{gene}.xlsx",
                              gene=gene,
                              threshold=threshold
                              )

            for gene in genes:
                for g1, g2 in combinations(group_keys, 2):
                    contig(df=df,
                           group1=g1,
                           group2=g2,
                           file_name=f'{folder_name}/{label_key}={label}--{g1}+{g2}-{gene}.xlsx',
                           gene=gene,
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
                    try:
                        sc.tl.rank_genes_groups(adata, vs, method='wilcoxon', use_raw=False)
                        dge = sc.get.rank_genes_groups_df(adata, 'positive', pval_cutoff=0.05).sort_values('logfoldchanges',
                                                                                                           ascending=False)
                        if dge.shape[0] < 2:
                            dge.to_excel(
                                f"{folder_name}/global-{label_key}={label}-{vs.split('_', 1)[0]}+_vs_{vs.split('_', 1)[0]}-NoGene.xlsx")
                        else:
                            dge.to_excel(
                                f"{folder_name}/global-{label_key}={label}-{vs.split('_', 1)[0]}+_vs_{vs.split('_', 1)[0]}-NoGene.xlsx")

                    except:
                        print(f"Couldn't calculated global DGE for {vs}.")
            cdata = adata.copy()
            for factor in factors:

                for lev in adata.obs[factor].unique():
                    adata = cdata.copy()
                    adata = adata[adata.obs[factor] == lev].copy()

                    for vs in versus:
                        try:
                            sc.tl.rank_genes_groups(adata, vs, method='wilcoxon', use_raw=False)
                            dge = sc.get.rank_genes_groups_df(adata, 'positive',
                                                              pval_cutoff=0.05).sort_values('logfoldchanges',
                                                                                           ascending=False)

                            if dge.shape[0] < 2:
                                dge.to_excel(
                                f"{folder_name}/{label_key}={label}+{factor}={lev}-{vs.split('_', 1)[0]}+_vs_{vs.split('_', 1)[0]}-_NoGene.xlsx")
                            else:
                                dge.to_excel(
                                f"{folder_name}/{label_key}={label}+{factor}={lev}-{vs.split('_', 1)[0]}+_vs_{vs.split('_', 1)[0]}-.xlsx")
                        except:
                            print(f"Couldn't calculated DGE for {vs} among {label} of {label_key} in group {lev} of {factor}")
    shutil.make_archive(folder_name, 'zip', folder_name)