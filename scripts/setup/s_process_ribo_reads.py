# If the file is too large to load in memory, you can split the sam file up and adjust the script
import sys
import h5py
import h5max
import pandas as pd
import numpy as np
from scipy import sparse
from tqdm import tqdm
import polars as pl

def process_ribo_reads(h5, ribo_metadata):
    ## Store Ribosome signal by experiment/read_length
    
    experiments = pd.read_csv(ribo_metadata, header=None, sep="\s+").values[:,0]
    header_dict = {2:'tr_ID', 3:'pos', 9:'read'}
    usecols = [2,3,9]

    f = h5py.File(h5, 'a')
    if 'riboseq' not in f['transcript'].keys():
        f['transcript'].create_group('riboseq')

    read_lens = np.arange(20,41)
    read_len_l = len(read_lens)

    tr_ids_f = np.array(f['transcript/id'])
    tr_lens_f = np.array(f['transcript/tr_len'])

    for experiment in experiments:
        print(f'Loading in {experiment}...')
        if experiment in f['transcript/riboseq'].keys():
            del f[f'transcript/riboseq/{experiment}']
        f['transcript/riboseq'].create_group(experiment)
        experiment_f = f[f'transcript/riboseq/{experiment}'].create_group('5')
        df = pl.read_csv(f'ribo/{experiment}/out/{experiment}_aligned_tran.sam', has_header=False, 
                        comment_char='@', columns=usecols, sep='\t')
        df.columns = list(header_dict.values())

        print('Filtering on read lens...')
        df = df.with_columns(
            pl.col("read").str.lengths().alias("read_len")
        )
        df = df.filter((pl.col('read_len') <= 40) & (pl.col('read_len') >= 20))
        mask_f = np.isin(tr_ids_f, df['tr_ID'].unique().to_numpy().astype('S'))

        print('Constructing empty datasets...')
        sparse_array = [sparse.csr_matrix((read_len_l, w)) for w in tqdm(tr_lens_f[~mask_f])]
        # construct riboseq data array
        riboseq_data = np.empty(len(mask_f), dtype=object)
        # fill with empty entries for which no data is available
        riboseq_data[~mask_f] = sparse_array
        df = df.sort('tr_ID')

        print('Mapping reads...')
        for tr_id, group in tqdm(df.groupby('tr_ID'), total=len(df['tr_ID'].unique())):
            mask_tr = tr_ids_f == tr_id.encode()
            tr_reads = np.zeros((read_len_l, tr_lens_f[mask_tr][0]), dtype=np.uint32)
            for row in group.rows():
                tr_reads[row[3]-20, row[1]-1] += 1
            riboseq_data[mask_tr] = sparse.csr_matrix(tr_reads)
        
        print('Saving data...')
        h5max.store_sparse(experiment_f, riboseq_data, format='csr')
        num_reads = [s.sum() for s in riboseq_data]
        experiment_f.create_dataset('num_reads', data=np.array(num_reads).astype(int))
    f.close()

if __name__ == "__main__":
    if len(sys.argv) > 3 or len(sys.argv) < 3:
        print("process_ribo_reads [h5_file] \n h5_file: path"
              " to hdf5 file containing transcriptome input data")
    else:
        print(sys.argv)
        process_ribo_reads(*sys.argv[1:])
