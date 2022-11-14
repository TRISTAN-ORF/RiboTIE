# Instructions

Following are the instructions on how to set up Riboformer and pre-process data. 
The requirements of the hardware are ca.

- \>100GB of RAM       (mostly data pre-processing)
- 500Gb of storage   (mostly dependent on RIBO data)
- 1 24Gb GPU
- 4 CPU's

## Setting up

For this project, a defined folder structure is used:

```
riboformer                              root folder
├── data                                input data
│   ├── genome                          reference genome
│   ├── ribo                            ribosome profiling
├── scripts                             execution scripts
├── models                              trained models
├── outputs                             model predictions

```
The process of setting up is mostly automated. Many of the following steps are accomponied by scripts. Download and unzip the folder containing scripts to a desired location, given access to the aforementioned hardware requirements.

#%LINK#

Note that scripts can take up to multiple hours, and require to be run from a terminal that will be active for that amount of time. `tmux` software can be used to detach from a terminal running a script without terminating it.

**In order to successfully perform the computations listed within the script files, it is important to run the scripts from within the `scripts/setup/` folder.**

### Software requirements

Several software packages need to be installed for the scripts to succeed. 

#### System

These packages need to be installed and accessible through your `PATH` variable:

- [gffreads](https://github.com/gpertea/gffread)
- [samtools](https://www.htslib.org/)
- [STAR](https://github.com/alexdobin/STAR)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [conda](https://docs.conda.io/en/latest/)

**Note**: It is possible to install `fastqc`, `cutadapt` and `gffreads` through `conda`, and include it as part of a python environment.

#### Python environment

The application of predictive tools and part of the data processing is performed through python. As with any python project, it is recommended to create a custom environment using `conda`. 

To create a new environment, run:

```
conda create --name riboformer_env
```

Activate the new environment to ensure new python packages are installed within:

```
conda activate riboformer_env
```

`PyTorch` is used as the deep learning library. Follow the instructions [here](https://pytorch.org/get-started/locally/) to install `PyTorch`. GPU support is highly recommended.



- gtfparse
- pandas
- tqdm
- polars
- pysam
- transcript_transformer


## Data Preprocessing

Data preprocessing includes setting up a reference genome and mapping ribosome profiling reads. 
After

### Setting up the reference genome

Start with setting up the reference genome. Running the script `scripts/setup/1_get_ref_genome.sh` will download and map the relevant genome into the `data/genome/` folder. Note that the version number of the human genome is a variable at the start of the script and can be altered in the future.
Within the script folder, run:

```
bash 1_get_ref_genome.sh
```

The next step is to create indexes for our mapping software `STAR`. Run:

```
bash 2_STAR_ref_genome.sh
```

### Mapping ribosome profiling data

To map the ribosome profiling data to the transcriptome, the ribosomal data need first be processed and mapped. Different ribosome experiments exist within the `data/ribo/` subfolder. For each experiment, create and name a folder after the experiment. Within the folder, place the `.fastq` ribosome file following the same name convention.

`data/ribo/metadata.txt` is a tab delimited metadata file containing the names (i.e. folder name) and adapter sequences for each experiment. Include the relevant information when adding new experiments in the folder.


An example folder layout:

```
riboformer                              root folder
├── data                                input data
│   ├── ribo                            ribosome profiling
│   │   ├── experiment_1                experiment folder
│   │   │   ├── experiment_1.fastq   
│   │   ├── experiment_2                experiment folder
│   │   │   ├── experiment_2.fastq
│   │   ├── metadata.txt                experiment metadata file

...
```

After adding all experiments and editing the metadata.txt file, run from within the `scripts/setup` folder:

```
bash 3_map_ribo_reads.sh
```

This will cut adapter sequences and map the reads to the genome and transcriptome.

**Note**: The script `3_map_ribo_reads.sh` will always process all experiments listed within the `metadata.txt` file. Consider using multiple `metadata.txt` files and editing the last line of the script (i.e. which metadata file is read) when it is not desired to start (re-)mapping the ribosome reads for all experiments.

### Formatting the model input data

The input data for the model is stored using the `hdf5` format. This file type allows quick access of selected data without the requirement of loading the full file into memory. It is furthermore faster than loading data from thousands of individual smaller files. the `hdf5` format supports storing data in a hierarchial data structure. 
Data saving and loading is achieved using the `h5py` python module.

The data is stored by transcript and information type. Information belonging to a single transcript is mapped according to the index they populate within each `h5py.dataset`.  Using this structure, it is fast to load various types of data given a transcript index. Ribosome reads are stored by read length and 5' position. Multiple experiments can be stored within a single `hdf5` file under the `/transcript/ribo/` file path. As ribosome mapping data is sparse, the python modules `scipy.sparse` and `h5max` are used to store and load sparse matrices (saved as multiple vector arrays under `data`, `indices`, `indptr` and `shape`). 

The internal file structure is as follows:

```
GRCh38_v107.h5                              (h5py.file)
    transcript                              (h5py.group)
    ├── tis                                 (h5py.dataset, dtype=vlen(int))
    ├── contig                              (h5py.dataset, dtype=str)
    ├── id                                  (h5py.dataset, dtype=str)
    ├── seq                                 (h5py.dataset, dtype=vlen(int))
    ├── ribo                                (h5py.group)
    │   ├── SRR0000001                      (h5py.group)
    │   │   ├── 5                           (h5py.group)
    │   │   │   ├── data                    (h5py.dataset, dtype=vlen(int))
    │   │   │   ├── indices                 (h5py.dataset, dtype=vlen(int))
    │   │   │   ├── indptr                  (h5py.dataset, dtype=vlen(int))
    │   │   │   ├── shape                   (h5py.dataset, dtype=vlen(int))
    │   ├── ...
    │   ....
    
```

In a first step, the `GRCh38_v107.h5` file is generated containing the transcriptome data (e.g. `tis`, `seq`, ...). This is achieved by parsing the reference `gtf` and `fasta` files from `data/genome/` to generate the `GRCh38_v107.h5` file within the `data/` folder. 

Run from the `script/setup/` folder:

```
bash 4_data_ref_genome.sh
```


The next step parses the mapped ribosome reads and incorporates them within the newly generated `GRCh38_v107.h5` file. Based on the size of the output `*.sam` file generated during the mapping step, this process can require >100Gb of RAM

```
bash 5_parse_ribo_reads.sh
```

**Note**: This script will process all ribosome experiments listed within the `data/ribo/metadata.txt` file. When it is desired to process only part of the experiments, edit either the python script or `metadata.txt` file.

## Model Training