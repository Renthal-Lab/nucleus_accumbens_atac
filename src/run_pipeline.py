import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
_stderr = sys.stderr
# null = open(os.devnull,'wb')
sys.version

# Working directory containing snATAC-seq and snRNA-seq data
work_dir = '/path/to/workdir/'
tmp_dir = '/path/to/tmp'


#make a directory for to store the processed scRNA-seq data.
if not os.path.exists(os.path.join(work_dir, 'scRNA')):
    os.makedirs(os.path.join(work_dir, 'scRNA'))


# Read in AnnData using scanpy
anndata_path = "/path/to/NAc_IdentsStored_MuSeurat.h5ad"
import scanpy as sc
adata = sc.read_h5ad(anndata_path)


adata


if not os.path.exists(os.path.join(work_dir, 'scATAC')):
    os.makedirs(os.path.join(work_dir, 'scATAC'))
    
import pandas as pd
nucl_acc_atac_metadata_path = '/n/scratch3/users/x/xl266/scenic_plus/nucleus_accumbens_Mar25/data/nucleus_accumbens_metadata.csv'
nucl_acc_metadata_anchored_atac = pd.read_csv(nucl_acc_atac_metadata_path, index_col=0)

cell_data = pd.DataFrame(nucl_acc_metadata_anchored_atac['Identities'])
cell_data['sample_id'] = 'nucl_acc'
cell_data

# Get chromosome sizes
import pyranges as pr
import requests
import pandas as pd

target_url='https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'

chromsizes=pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
chromsizes


# Run CisTopic
bc_passing_filters = {'nucl_acc':[]}
bc_passing_filters['nucl_acc'] = list(set(cell_data.index))

# Load ATAC data
import pickle
fragments_path = '/n/scratch3/users/p/pab1164/NAc/fragments.tsv.gz'
fragments_dict = {'nucl_acc': fragments_path}
path_to_regions = {'nucl_acc':os.path.join(work_dir, 'scATAC/nucl_acc_peaks_atac.bed')}
path_to_blacklist= os.path.join(work_dir, 'data/mm10-blacklist.v2.bed')

# Create CisTopic object
from pycisTopic.cistopic_class import *
key = 'nucl_acc'
cistopic_obj = create_cistopic_object_from_fragments(
                            path_to_fragments=fragments_dict[key],
                            path_to_regions=path_to_regions[key],
                            path_to_blacklist=path_to_blacklist,
                            valid_bc=list(set(bc_passing_filters[key])),  # deleted ' & set(scRNA_bc)', 'metrics = '
                            n_cpu=1,
                            project=key,
                            split_pattern='-')
cistopic_obj.add_cell_data(cell_data, split_pattern='-')
print(cistopic_obj)
# Save CisTopic object
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/nucl_acc_cistopic_obj.pkl'), 'wb'))


# Topic modeling
import pickle
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/nucl_acc_cistopic_obj.pkl'), 'rb'))
from pycisTopic.cistopic_class import *
models=run_cgs_models(cistopic_obj,
                    n_topics=[2,4,10,16,32,48],
                    n_cpu=5,
                    n_iter=500,
                    random_state=555,
                    alpha=50,
                    alpha_by_topic=True,
                    eta=0.1,
                    eta_by_topic=False,
                    save_path=None,
                    _temp_dir = os.path.join(tmp_dir + 'ray_spill'))


# Save topic modeling results

if not os.path.exists(os.path.join(work_dir, 'scATAC/models')):
    os.makedirs(os.path.join(work_dir, 'scATAC/models'))

pickle.dump(models,
            open(os.path.join(work_dir, 'scATAC/models/nucl_acc_models_500_iter_LDA.pkl'), 'wb'))



# Evaluate models
models = pickle.load(open(os.path.join(work_dir, 'scATAC/models/nucl_acc_models_500_iter_LDA.pkl'), 'rb'))
cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/nucl_acc_cistopic_obj.pkl'), 'rb'))
from pycisTopic.lda_models import *

# Let pycisTopic pick the optimal number of topics. In this case, 32
model = evaluate_models(models,
                      select_model=None,
                      return_model=True,
                      metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                      plot_metrics=False,
                      save = os.path.join(work_dir, 'scATAC/models/evaluation.pdf'))


# Add chosen model to cistopic object
cistopic_obj.add_LDA_model(model)
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/nucl_acc_cistopic_obj.pkl'), 'wb'))


# Visualize topics
from pycisTopic.clust_vis import *
run_umap(cistopic_obj, target  = 'cell', scale=True)
plot_metadata(cistopic_obj, reduction_name = 'UMAP', variables = ['Identities'], save = os.path.join(work_dir, 'scATAC/topics_umap.pdf'))


# plot cell-topic probabilities
plot_topic(cistopic_obj, reduction_name = 'UMAP', num_columns = 4, save = os.path.join(work_dir, 'scATAC/topics_cell_topic_prob_umap.pdf'))


# Infer candidate enhanceres
import pickle
from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

# Calculate DARs per celltype
from pycisTopic.diff_features import *
imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable='Identities', var_features=variable_regions, split_pattern = '-')

# Save results
if not os.path.exists(os.path.join(work_dir, 'scATAC/candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, 'scATAC/candidate_enhancers'))
import pickle
pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'wb'))


Motif enrichment analysis using pycistarget

Load candidate enhancers
import pickle
region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
markers_dict = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/markers_dict.pkl'), 'rb'))


# Convert to dictionary of pyranges objects
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))
    
    
for key in region_sets.keys():
    print(f'{key}: {region_sets[key].keys()}')
    

motif_annotation = '/n/scratch3/users/x/xl266/scenic_plus/drg_atac_rna_v2/data/motif_data_scenic/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl'
encode_db_scores_path = '/n/scratch3/users/x/xl266/scenic_plus/drg_atac_rna_v2/data/motif_data_scenic/mm10_screen_v10_clust.regions_vs_motifs.scores.feather'
encode_db_rankings_path = '/n/scratch3/users/x/xl266/scenic_plus/drg_atac_rna_v2/data/motif_data_scenic/mm10_screen_v10_clust.regions_vs_motifs.rankings.feather'

if not os.path.exists(os.path.join(work_dir, 'motifs')):
    os.makedirs(os.path.join(work_dir, 'motifs'))

from scenicplus.wrappers.run_pycistarget import run_pycistarget
run_pycistarget(
    region_sets = region_sets,
    species = 'mus_musculus',
    save_path = os.path.join(work_dir, 'motifs'),
    ctx_db_path = encode_db_rankings_path,  # using renamed version (columns are region strings)
    dem_db_path = encode_db_scores_path,    # using renamed version (columns are region strings)
    annotation = ['Direct_annot', 'Orthology_annot', 'Motif_similarity_annot', 'Motif_similarity_and_Orthology_annot'],    # Added Mar 19, 2023
    path_to_motif_annotations = motif_annotation,
    run_without_promoters = True,
    n_cpu = 12,
    _temp_dir = os.path.join(tmp_dir, 'ray_spill'),
    annotation_version = 'v10nr_clust'
    )



# inferring enhancer-driven Gene Regulatory Networks (eGRNs) using SCENIC+
import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
import pickle

adata = sc.read_h5ad(anndata_path)    
cistopic_obj = dill.load(open(os.path.join(work_dir, 'scATAC/nucl_acc_cistopic_obj.pkl'), 'rb'))
menr = dill.load(open(os.path.join(work_dir, 'motifs/menr.pkl'), 'rb'))

# Mar 22, 2023: the .X slot of reclustered AnnData is actually normalized, so we replace w/ raw counts sparse matrix
adata.X = adata.layers['counts'].copy()
adata.raw = adata

# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')


from scenicplus.scenicplus_class import create_SCENICPLUS_object
import numpy as np
import scipy as sp

scplus_obj = create_SCENICPLUS_object(
    GEX_anndata = adata.raw.to_adata(),
    cisTopic_obj = cistopic_obj,
    menr = menr,
    multi_ome_mode = False,
    key_to_group_by = 'Identities',
    meta_cell_split = '.'     # not '_' since that appears in cell type labels
)

scplus_obj.X_EXP = sp.sparse.csr_matrix(scplus_obj.X_EXP)
scplus_obj.X_EXP = np.array(scplus_obj.X_EXP.todense())


# check with which biomart host our gene names match best

ensembl_version_dict = {'105': 'http://www.ensembl.org',
                        '104': 'http://may2021.archive.ensembl.org/',
                        '103': 'http://feb2021.archive.ensembl.org/',
                        '102': 'http://nov2020.archive.ensembl.org/',
                        '101': 'http://aug2020.archive.ensembl.org/',
                        '100': 'http://apr2020.archive.ensembl.org/',
                        '99': 'http://jan2020.archive.ensembl.org/',
                        '98': 'http://sep2019.archive.ensembl.org/',
                        '97': 'http://jul2019.archive.ensembl.org/',
                        '96': 'http://apr2019.archive.ensembl.org/',
                        '95': 'http://jan2019.archive.ensembl.org/',
                        '94': 'http://oct2018.archive.ensembl.org/',
                        '93': 'http://jul2018.archive.ensembl.org/',
                        '92': 'http://apr2018.archive.ensembl.org/',
                        '91': 'http://dec2017.archive.ensembl.org/',
                        '90': 'http://aug2017.archive.ensembl.org/',
                        '89': 'http://may2017.archive.ensembl.org/',
                        '88': 'http://mar2017.archive.ensembl.org/',
                        '87': 'http://dec2016.archive.ensembl.org/',
                        '86': 'http://oct2016.archive.ensembl.org/',
                        '80': 'http://may2015.archive.ensembl.org/',
                        '77': 'http://oct2014.archive.ensembl.org/',
                        '75': 'http://feb2014.archive.ensembl.org/',
                        '54': 'http://may2009.archive.ensembl.org/'}

import pybiomart as pbm
def test_ensembl_host(scplus_obj, host, species):
    dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Chromosome'] = annot['Chromosome'].astype('str')
    filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    gene_names_release = set(annot['Gene'].tolist())
    ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
    print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
    return ov

n_overlap = {}
for version in ensembl_version_dict.keys():
    print(f'host: {version}')
    try:
        n_overlap[version] =  test_ensembl_host(scplus_obj, ensembl_version_dict[version], 'mmusculus')
    except:
        print('Host not reachable')
v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")


biomart_host = str(ensembl_version_dict[v])


from scenicplus.wrappers.run_scenicplus import run_scenicplus

# TF file from https://resources.aertslab.org/cistarget/tf_lists/allTFs_mm.txt
# import ray
# ray.shutdown()

# Only use save_partial = True if you are running from the very beginning / earlier in the pipeline.
# Otherwise, scenicplus will resave the entire object (takes time) at every if save_partial statement

bedToBigBed_dir = '/n/scratch3/users/x/xl266/scenic_plus/drg_atac_rna_v2/'
         
try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['Identities'],
        species = 'mmusculus',
        assembly = 'mm10',
        tf_file = '/n/scratch3/users/x/xl266/scenic_plus/drg_atac_rna/data/allTFs_mm.txt',
        save_path = os.path.join(work_dir, 'scenicplus'),
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = True,
        export_to_loom_file = False,
        export_to_UCSC_file = True,
        path_bedToBigBed = bedToBigBed_dir,
        n_cpu = 20,
        _temp_dir = os.path.join(tmp_dir, 'ray_spill'))
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(work_dir, 'scenicplus/scplus_obj.pkl'), 'wb'), protocol=-1)
    raise(e)


