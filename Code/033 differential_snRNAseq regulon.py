#!/usr/bin/env python
# coding: utf-8


import scanpy as sc
import os
import glob
import pickle
import pandas as pd
import numpy as np
from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

DATA_FOLDER="/ChP_RNA/Scenic/results"
RESOURCES_FOLDER="/ChP_RNA/Scenic/resources"
DATABASE_FOLDER = "/ChP_RNA/Scenic/databases/"
SCHEDULER="123.122.8.24:8786"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "*.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
HS_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'allTFs_hg38.txt')
REGULONS_FNAME = os.path.join(DATA_FOLDER, "regulons.p")
MOTIFS_FNAME = os.path.join(DATA_FOLDER, "motifs.csv")



#load data
adata = sc.read_h5ad("/ChP_RNA/Scenic/results/Ex_ChP_sub.h5ad")

ex_matrix = adata.X.tocsc() if hasattr(adata.X, 'tocsc') else csc_matrix(adata.X)
ex_matrix = ex_matrix.T 

gene_names = adata.var_names.tolist()
tf_names = load_tf_names(HS_TFS_FNAME)
print(ex_matrix.shape)
print(len(gene_names))
print(len(tf_names))

row_attrs = {"Gene": np.array(gene_names)}
col_attrs = {"CellID": np.array(adata.obs_names)}

lp.create("/ChP_RNA/Scenic/results/Ex_ChP_sub.loom", ex_matrix, row_attrs, col_attrs)




#Set up dask
#from dask.distributed import Client, LocalCluster
#try:
#    client.close()
#except:
#    pass
    
#cluster = LocalCluster(n_workers=16, threads_per_worker=8, memory_limit='64GB')
#client = Client(cluster)

#Dask is terrible at processing anything in parralel unless your set is the very specific use case it works with. Seriously stop recommending it as default.
#adjacencies = grnboost2(expression_data=ex_matrix, gene_names=gene_names, tf_names=tf_names, verbose=True)

get_ipython().system('arboreto_with_multiprocessing.py      /ChP_RNA/Scenic/results/Ex_ChP_sub.loom      /ChP_RNA/Scenic/resources/allTFs_hg38.txt      --method grnboost2      --output /ChP_RNA/Scenic/results/adjacencies.tsv      --num_workers 32      --seed 777')




adata = sc.read_h5ad("/ChP_RNA/Scenic/results/Ex_ChP_full.h5ad")
row_attrs = {"Gene": np.array(adata.var_names) ,}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}

lp.create("/ChP_RNA/Scenic/results/Ex_ChP.loom", adata.X.transpose(), row_attrs, col_attrs)




import glob
# ranking databases
f_db_glob = "/ChP_RNA/Scenic/databases/*feather"
f_db_names = ' '.join( glob.glob(f_db_glob) )

# motif databases
f_motif_path = "/ChP_RNA/Scenic/resources/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

get_ipython().system('pyscenic ctx /ChP_RNA/Scenic/results/adjacencies.tsv      {f_db_names}      --annotations_fname {f_motif_path}      --expression_mtx_fname /ChP_RNA/Scenic/results/Ex_ChP.loom      --output /ChP_RNA/Scenic/results/reg.csv      --nes_threshold 2      --rank_threshold 10000      #--mode custom_multiprocessing      --min_genes 5      --num_workers 40')



import pandas as pd
import loompy as lp

adjacencies = pd.read_csv("/ChP_RNA/Scenic/results/adjacencies.tsv", sep='\t')

with lp.connect("/ChP_RNA/Scenic/results/Ex_ChP.loom", mode='r') as lf:
    ex_matrix = pd.DataFrame(lf[:, :], index=lf.ra.Gene, columns=lf.ca.CellID).T
    
modules = list(modules_from_adjacencies(adjacencies, ex_matrix, rho_threshold=0.01, rho_mask_dropouts=False)) #!!!!!Setting to F may cause errors

db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
len(dbs)




with ProgressBar():
   df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, client_or_address = 'custom_multiprocessing', num_workers=16)

# Create regulons from this table of enriched motifs.
regulons = df2regulons(df)

# Save the enriched motifs and the discovered regulons to disk.
df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)




print(f"Total regulons: {len(regulons)}")
print(f"First regulon object: {regulons[0]}")
for reg in regulons[:5]:
    print(f"TF: {reg.name:15} | Targets: {len(reg.genes)}")




auc_mtx = aucell(ex_matrix, regulons, num_workers=16)
auc_mtx.to_csv("/ChP_RNA/Scenic/results/auc_mtx_.csv")
auc_mtx = auc_mtx.T #transpose, genes are often cols rather than rows
auc_mtx.to_csv("/ChP_RNA/Scenic/results/auc_mtx_Transposed.csv")




import pandas as pd
import numpy as np
from pyscenic.binarization import binarize

auc_mtx = pd.read_csv("/ChP_RNA/Scenic/results/auc_mtx_.csv", index_col=0)

if auc_mtx.iloc[:, 0].dtype == 'object':
    auc_mtx = auc_mtx.set_index(auc_mtx.columns[0])

auc_mtx = auc_mtx.apply(pd.to_numeric)

binary_mtx, thresholds = binarize(auc_mtx)
#relax thresholds
scale = 0.7
thresholds_relaxed = {k: v * scale for k, v in thresholds.items()}
binary_mtx = auc_mtx.gt(pd.Series(thresholds_relaxed))

binary_mtx.to_csv("/ChP_RNA/Scenic/results/binary_regulon_activity.csv")



import pickle

with open("/ChP_RNA/Scenic/results/regulons.p", "rb") as f:
    regulons = pickle.load(f)

foxj1 = [r for r in regulons if r.name.startswith("FOXJ1")][0]

print(foxj1)




import pandas as pd

edges = []

for reg in regulons:

    tf = reg.transcription_factor

    for gene, weight in reg.gene2weight.items():

        edges.append({
            "TF": tf,
            "Target": gene,
            "Weight": weight,
            "NES": reg.nes,
            "Regulon": reg.name
        })

edges = pd.DataFrame(edges)

edges.to_csv("/ChP_RNA/Scenic/results/SCENIC_regulon_edges.csv", index=False)

