import pandas as pd
import scanpy as sc
from anndata import AnnData
from pandas import DataFrame
from itertools import combinations
from collections import OrderedDict
import os
import shutil


def _check_expression_mtx(adata):
    m = adata.X < 0
    if m.size > 0:
        return False


def _check_folder(folder_name):
    if not os.path.isdir(folder_name):
        os.mkdir(folder_name)


def get_annotated(adata, genes):
    if _check_expression_mtx(adata):
        print('Expression matrix seems to be Z transformed!!!')

    if isinstance(genes, list):
        # Add expression levels to df
        subdata = adata[:, adata.var_names.isin(genes)].copy()
        df = subdata.to_df()
        df = pd.concat([subdata.obs, df], axis=1)

        # Annotate Positivity
        for gene in genes:
            gname = f'{gene}_status'
            df[gname] = df[gene].apply(lambda x: 'positive' if x > 0 else 'negative')

        # Add annotation to adata
        for gene in genes:
            gname = f'{gene}_status'
            adata.obs[gname] = df[gname].copy()

        return df, adata


def pct_table(df,
              group,
              gene,
              file_name='table.xlsx',
              threshold=0  # For this to be valid data shouldn't be zero centered.
              ):
    vals = OrderedDict({})
    for cell in df[group].value_counts().sort_index().index:
        allcells = df[df[group] == cell].shape[0]
        cds = df[(df[group] == cell) & (df[gene] > threshold)].shape[0]
        if allcells == 0:
            pct = 0
        else:
            pct = round(cds / allcells * 100, 2)
        vals[cell] = [allcells, cds, pct]
    table = pd.DataFrame(vals, index=['All Cells', f'{gene} Expressing', 'Percentage'])
    table.to_excel(file_name)
    return table


def contig(df: DataFrame,
           group1,
           group2,
           gene,
           file_name='table.xlsx',
           threshold=0
           ):
    grp = OrderedDict({})
    for i in df[group1].value_counts().sort_index().index:
        grp[i] = pct_table(df[df[group1] == i], group2, file_name=file_name, gene=gene, threshold=threshold)
    table = pd.concat(grp)
    table.to_excel(file_name)
    return table


def analyze_pct(data: DataFrame,
                label_keys: list,
                factor_keys: list,
                genes: list,
                analyze_global: bool = True,
                threshold=0,
                folder_name: str = 'pct'
                ):
    _check_folder(folder_name)

    for label_key in label_keys:

        if analyze_global:
            for gen in genes:
                pct_table(data,
                          group=label_key,
                          gene=gen,
                          file_name=f"{folder_name}/global-{label_key}-{gen}.xlsx")
        else:
            pass

        dff = data.copy()
        for gen in genes:
            for factor_key in factor_keys:
                contig(df=dff,
                       group1=label_key,
                       group2=factor_key,
                       file_name=f'{folder_name}/contingency--{label_key}+{factor_key}-{gen}.xlsx',
                       gene=gen,
                       threshold=threshold)

            if len(factor_keys) > 1:
                _check_folder(f'{folder_name}/detailed/')
                for lbs in label_keys:
                    _check_folder(f'{folder_name}/detailed/{lbs}')

                for factor1, factor2 in combinations(factor_keys, 2):

                    for label_level in dff[label_key].unique():
                        celldf = dff.copy()
                        celldf = celldf[celldf[label_key] == label_level].copy()
                        contig(df=celldf,
                               group1=factor1,
                               group2=factor2,
                               file_name=f'{folder_name}/detailed/{label_key}/contingency--{label_key}={label_level}\
                               --{factor1}+{factor2}-{gen}.xlsx',
                               gene=gen,
                               threshold=threshold
                               )

    shutil.make_archive(folder_name, 'zip', folder_name)


def analyze_dge(adata: AnnData,
                label_keys: list,
                factor_keys: list,
                compare: list,
                analyze_global: bool = True,
                analyze_interaction: bool = True,
                folder_name: str = 'dge'):
    _check_folder(folder_name)

    # Iterate over labels
    for label_key in label_keys:
        for label in adata.obs[label_key].unique():
            tempdata_label = adata[adata.obs[label_key] == label].copy()

            if analyze_global:
                for vs in compare:
                    try:
                        sc.tl.rank_genes_groups(tempdata_label, vs, method='wilcoxon', use_raw=False)
                        dge = sc.get.rank_genes_groups_df(tempdata_label, 'positive', pval_cutoff=0.05).sort_values(
                            'logfoldchanges',
                            ascending=False)
                        if dge.shape[0] < 2:
                            dge.to_excel(
                                f"{folder_name}/global-{label_key}={label}-{vs.split('_', 1)[0]}+\
                                _vs_{vs.split('_', 1)[0]}-NoGene.xlsx")
                        else:
                            dge.to_excel(
                                f"{folder_name}/global-{label_key}={label}-{vs.split('_', 1)[0]}+\
                                _vs_{vs.split('_', 1)[0]}-.xlsx")

                    except:
                        print(f"Couldn't calculated global DGE for {vs} in {label}.")

            for factor_key in factor_keys:
                for lev in adata.obs[factor_key].unique():
                    tempdata_level = tempdata_label[tempdata_label.obs[factor_key] == lev].copy()

                    for vs in compare:
                        try:
                            sc.tl.rank_genes_groups(tempdata_level, vs, method='wilcoxon', use_raw=False)
                            dge = sc.get.rank_genes_groups_df(tempdata_level, 'positive',
                                                              pval_cutoff=0.05).sort_values('logfoldchanges',
                                                                                            ascending=False)

                            if dge.shape[0] < 2:
                                dge.to_excel(
                                    f"{folder_name}/{label_key}={label}+{factor_key}={lev}\
                                    -{vs.split('_', 1)[0]}+_vs_{vs.split('_', 1)[0]}-_NoGene.xlsx")
                            else:
                                dge.to_excel(
                                    f"{folder_name}/{label_key}={label}+{factor_key}={lev}-{vs.split('_', 1)[0]}+\
                                    _vs_{vs.split('_', 1)[0]}-.xlsx")
                        except:
                            print(
                                f"Couldn't calculated DGE for {vs} in the group {label}--{factor_key}={lev}")

                    if len(factor_keys) > 1 and analyze_interaction:
                        other_factors = [x for x in factor_keys if x != factor_key]
                        for other_factor in other_factors:
                            for other_level in tempdata_level.obs[other_factor].unique():
                                tempdata_level_interact = tempdata_level[
                                                          tempdata_level.obs[other_factor] == other_level, :].copy()

                                for vs in compare:
                                    try:
                                        sc.tl.rank_genes_groups(tempdata_level_interact, vs, method='wilcoxon',
                                                                use_raw=False)
                                        dge = sc.get.rank_genes_groups_df(tempdata_level_interact, 'positive',
                                                                          pval_cutoff=0.05).sort_values(
                                            'logfoldchanges',
                                            ascending=False)

                                        if dge.shape[0] < 2:
                                            dge.to_excel(
                                                f"{folder_name}/{label_key}={label}+{factor_key}={lev}+{other_factor}\
                                                ={other_level}-{vs.split('_', 1)[0]}+_vs_{vs.split('_', 1)[0]}-_NoGene.xlsx")
                                        else:
                                            dge.to_excel(
                                                f"{folder_name}/{label_key}={label}+{factor_key}={lev}+{other_factor}\
                                                ={other_level}-{vs.split('_', 1)[0]}+_vs_{vs.split('_', 1)[0]}-.xlsx")
                                    except:
                                        print(
                                            f"Couldn't calculated DGE for {vs} in the group : {label}--{factor_key}\
                                            ={lev}--{other_factor}={other_level}")

    shutil.make_archive(folder_name, 'zip', folder_name)
