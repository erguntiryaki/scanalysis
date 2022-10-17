import scipy as ir
import matplotlib.pyplot as plt
from anndata import AnnData

def donut(df, group, title='', save=None):
    lbs = df[group].value_counts().sort_index().index.to_list()
    sz = df[group].value_counts().sort_index().values
    try:
        explode = [0.05 for i in range(df[group].nunique())]
        _ = explode.pop(0)
        explode.insert(0, 0.)

        plt.pie(sz, labels=lbs, autopct='%1.2f%%', pctdistance=0.78,
                explode=explode)
    except ValueError:
        plt.pie(sz, labels=lbs, autopct='%1.2f%%', pctdistance=0.78)
    centre_circle = plt.Circle((0, 0), 0.55, fc='white')
    fig = plt.gcf()

    fig.gca().add_artist(centre_circle)
    plt.title(title)
    plt.tight_layout()
    if save is not None:
        plt.savefig(f'{save}.png', dpi=300)

def filter_vdjdb(vdjdb=None, species='HomoSapiens', mhc_class='MHCI', min_length=6):
    if vdjdb is not None:
        ir.dataets.vdjdb()
    vdjdb = vdjdb[(vdjdb.obs['species'] == species) &
                  (vdjdb.obs['mhc.class'] == mhc_class) &
                  (~ vdjdb.obs.IR_VDJ_1_junction_aa.isna()), :]

    arr = vdjdb.obs['IR_VDJ_1_junction_aa'].apply(lambda x: False if len(str(x)) < min_length else True)
    vdjdb = vdjdb[arr.to_list(), :].copy()
    return vdjdb

class Experiment:
    def __init__(self, adata, gene):
        self.adata = adata,
        self.gene = gene,
        self.adata_pos = adata[adata.obs.CD20_status == 'positive']
        self.adata_neg = adata[adata.obs.CD20_status == 'negative']



