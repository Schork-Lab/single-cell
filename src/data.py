import os

import pandas as pd

from helpers import *


def load_mapped_data(idx_dir=None):
    if not idx_dir:
        idx_dir = os.path.join(local_data, 'nov_24_run',
                             'HBPC_RSEM_RNASeq_PV_12012014',
                             'BAM_files', 'genome')

    mapped_stats = {}
    for filename in os.listdir(idx_dir):
        if filename.endswith('.idxstats'):
            sample = filename.split('.')[0]
            filename = os.path.join(idx_dir, filename)
            df = pd.read_table(filename, sep='\t',
                               names=['contig', 'length', 'reads_mapped', 'reads_unmapped'])
            ercc_mapped = df[df.contig.str.contains('ERCC')]['reads_mapped'].sum()
            other_mapped = df[~df.contig.str.contains('ERCC')]['reads_mapped'].sum()
            unmapped = df['reads_unmapped'].sum()
            mapped_stats[sample] = (ercc_mapped, other_mapped, unmapped)
    mapped_df = pd.DataFrame(mapped_stats).T
    mapped_df.columns = ['ERCC', 'GENOME', 'UNMAPPED']
    return mapped_df

def load_tpms(data_dir=None, sample_map=None):
    '''
    
    :param data_dir:
    :type data_dir:
    :param sample_map:
    :type sample_map:
    '''
    if not data_dir:
        data_dir = os.path.join(local_data, 'nov_24_run', 'HBPC_RSEM_RNASeq_PV_12012014',
                                                        'expression_data', 'genes')
    if not sample_map:
        sample_map = get_sample_map()

    sample_dfs = []
    for fn in os.listdir(data_dir):
        sample_id = fn.split('.')[0]
        if fn.endswith('.genes.results') and sample_id in sample_map:
            sample_name = sample_map[sample_id]
            sample_df = read_sample(os.path.join(data_dir, fn),
                                                                      sample_name)
            sample_dfs.append(sample_df)

    tpm_df = pd.concat(sample_dfs, axis=1)

    return tpm_df

def paperify_df(df, sample_map=None):
    '''
    Converts the columns of a df to match the names that should be on the paper.
    '''
    if not sample_map:
        sample_map = get_sample_map()
    df = df[[column for column in df.columns
                    if column in sample_map ]]
    df.columns = [sample_map[column] for column in df.columns]
    return df


def read_sample(filename, sample_name=None):
    '''
    
    :param filename:
    :type filename:
    :param sample_name:
    :type sample_name:
    '''

    if not sample_name:
        sample_name = os.path.basename(filename).split('.')[0]
    df = pd.read_table(filename, index_col=0)
    df[sample_name] = df['TPM']
    return df[[sample_name]]

def filter_df(tpm_df, genes_only=True,
                         neuronal=True, non_neuronal=True, pooled=True):
    '''
    
    :param tpm_df:
    :type tpm_df:
    :param genes_only:
    :type genes_only:
    :param neuronal:
    :type neuronal:
    :param non_neuronal:
    :type non_neuronal:
    :param pooled:
    :type pooled:
    '''

    if genes_only:
        tpm_df = tpm_df.ix[[index for index in tpm_df.index
                                                     if index.find('ENSG') != -1]]

    columns = []
    if neuronal:
        columns += [column for column in tpm_df.columns
                                     if column.find('Neuronal nuclei') != -1]
    if non_neuronal:
        columns += [column for column in tpm_df.columns
                                     if column.find('Non-neuronal nuclei') != -1]
    if pooled:
        columns += [column for column in tpm_df.columns
                                     if column.find('Total RNA') != -1]

    return tpm_df[columns].sort(axis=1)

def get_low_mid_high_genes(control_tpm, cutoffs=(0.33, 0.66)):
    '''
    
    :param control_tpm:
    :type control_tpm:
    :param cutoffs:
    :type cutoffs:
    '''
    control_expressed = control_tpm[control_tpm > 0]
    low_cutoff = control_expressed <= control_expressed.quantile(cutoffs[0])
    mid_cutoff = ((control_expressed > control_expressed.quantile(cutoffs[0])) &
                                  (control_expressed < control_expressed.quantile(cutoffs[1])))
    high_cutoff = (control_expressed >= control_expressed.quantile(cutoffs[1]))
    low_expressed = control_expressed[low_cutoff]
    mid_expressed = control_expressed[ mid_cutoff]
    high_expressed = control_expressed[high_cutoff]
    return low_expressed, mid_expressed, high_expressed
