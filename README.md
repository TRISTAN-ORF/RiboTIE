<div align="center">
<h1>🧮 RIBO-former</h1>

*Driving coding sequence discovery since 2023*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8059764.svg)](https://doi.org/10.5281/zenodo.8059764)
[![GitHub license](https://img.shields.io/github/license/jdcla/RIBO_former)](https://github.com/jdcla/RIBO_former/blob/main/LICENSE.md)
[![GitHub issues](https://img.shields.io/github/issues/jdcla/RIBO_former)](https://github.com/jdcla/RIBO_former/issues)
[![GitHub stars](https://img.shields.io/github/stars/jdcla/RIBO_former)](https://github.com/jdcla/RIBO_former/stargazers)



</div>

## 📋 About
[RIBO-former](https://doi.org/10.1101/2023.06.20.545724) is created to annotate translation initiation sites on transcripts using ribosome profiling data. This repository contains the instructions to run RIBO-former on custom data.

The data, model parameters, and benchmark data from **the article** are featured in a [separate repository](https://github.com/jdcla/RIBO_former_paper).

When interested in more advanced features, such as using a custom transformer architecture, we refer the user manual of the [transcript-transformer package](https://github.com/jdcla/transcript_transformer), created in support of this tool. 

Make sure to check out [TIS Transformer](https://github.com/jdcla/TIS_transformer) as well, a similar tool for the delineation of novel coding sequences using transcript sequence data rather than ribosome profiling data.

## 📖 User guide

Following are the instructions on how to set up RIBO-former and pre-process data. 
The hardware requirements are:

- ~100-150GB of RAM       (mostly data pre-processing)
- ~500GB of storage       (mostly dependent on RIBO data)
- 1 24GB vRAM GPU
- ~4 CPU's

For now, a defined folder structure is used:

```
RIBOformer_tool                         root folder
├── data                                input data
│   ├── genome                          reference genome
│   ├── ribo                            ribosome profiling
├── scripts                             execution scripts
│   ├── setup                           setup scripts
│   ├── train                           train scripts
├── models                              trained models
├── outputs                             model predictions

```
The process of setting up is mostly automated. Many of the following steps are accomponied by scripts. Download and unzip the folder containing scripts to a desired location, given access to the aforementioned hardware requirements.

The folder and scripts can be downloaded using `git clone`

Note that scripts can take up to multiple hours, and require to be run from a terminal that will be active for that amount of time. `tmux` software can be used to detach from a terminal running a script without terminating it.

**Note**: In order to successfully perform the computations listed within the script files, it is important to run the scripts from within the `scripts/setup/` folder. Work is being done to integrate all steps into a single end-to-end pipeline**

### Software requirements

#### System
Several software packages need to be installed for the scripts to succeed. 

These packages need to be installed and accessible through your `PATH` variable:

- [samtools](https://www.htslib.org/)
- [STAR](https://github.com/alexdobin/STAR)
- [cutadapt](https://cutadapt.readthedocs.io/en/stable/)
- [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [pyfaidx](https://pypi.org/project/pyfaidx/#description)

**Note**: It is possible to install `fastqc`, `cutadapt` and `pyfaidx` through `conda`/`pip`, and include it as part of a python environment.

#### Python environment

The application of predictive tools and part of the data processing is performed through python. As with any python project, it is recommended to create a custom environment (e.g. using `conda`). Package installation is often easier achieved using `pip`, as conda often finds conflicts between packages that are typically not problematic.

To create a new environment, run:

```
conda create --name riboformer_env
```

Activate the new environment to ensure new python packages are installed within:

```
conda activate riboformer_env
```

`PyTorch` is used as the deep learning library. Follow the instructions [here](https://pytorch.org/get-started/locally/) to install `PyTorch` first. GPU support is necessary. 

The following python packages are used:
- gtfparse
- pyfaidx
- pandas
- tqdm
- polars
- transcript_transformer

Install using pip (ensure you have activated the conda environment):

```bash
pip install gtfparse pandas tqdm polars pysam transcript_transformer
```

Make sure all the steps are run with the custom evironment activated.

### Data Preprocessing

Data preprocessing includes setting up a reference genome and mapping ribosome profiling reads. 
After

#### Setting up the reference genome

Start with setting up the reference genome. Running the script `scripts/setup/1_get_ref_genome.sh` will download and map the relevant genome into the `data/genome/` folder. Note that the version number of the human genome is a variable at the start of the script and can be altered in the future.
Within the script folder, run:

```
bash 1_get_ref_genome.sh
```

The next step is to create indexes for our mapping software `STAR`. Run:

```
bash 2_STAR_ref_genome.sh
```

#### Mapping ribosome profiling data

To map the ribosome profiling data to the transcriptome, the ribosomal data need first be processed and mapped. Different ribosome experiments exist within the `data/ribo/` subfolder. For each experiment, create and name a folder after the experiment. Within the folder, place the `.fastq` ribosome file following the same name convention.

`data/ribo/metadata.txt` is a tab delimited metadata file containing the names (i.e. folder name) and adapter sequences for each experiment. Include the relevant information when adding new experiments in the folder.


An example folder layout:

```
RIBOformer_tool                         root folder
├── data                                input data
│   ├── ribo                            ribosome profiling
│   │   ├── experiment_1                experiment folder
│   │   │   ├── experiment_1.fastq   
│   │   ├── experiment_2                experiment folder
│   │   │   ├── experiment_2.fastq
│   │   ├── metadata.txt                experiment metadata file

...
```

After adding all experiments and editing the `metadata.txt` file, run from within the `scripts/setup` folder:

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
bash 4_process_transcriptome.sh
```


The next step parses the mapped ribosome reads and incorporates them within the newly generated `GRCh38_v107.h5` file. Based on the size of the output `*.sam` file generated during the mapping step, this process can require >100Gb of RAM

```
bash 5_process_ribo.sh
```

**Note**: This script will process all ribosome experiments listed within the `data/ribo/metadata.txt` file. When it is desired to process only part of the experiments, edit either the python script or `metadata.txt` file.

### Model Training

The scripts used to train the model are located under `/scripts/train`. The [transcript_transformer](https://github.com/jdcla/transcript_transformer) python package handles training and use of the predictive models. 

RIBO-former can be used to map the full translatome of a ribosome profiling experiment. For this, multiple models are used that are trained on different folds of the data. The parts of the transcriptome excluded from the training and model selection process can be used for mapping (i.e. unbiased). In the following set-up, we use 5/6 of the data to train/select a model with which the remaining 1/6 of the transcriptome is used to map the TIS on.

The script `/scripts/train/train.sh` executes the training of six models trained on the various folds. Each model takes around ~12 hours to converge. It is possible to run the script as is, training all models in series over a total time of around 3 days. Alternatively, if multiple GPUs are available, it might be preferable to break the script up and train the models in parallel.

When `transcript_transformer` is called, it requires a dictionary file as input that lays out the input file structure of the `hdf5` file and the data used to train the model. This file is located under `scripts/train/template.json`. **Note that it is important to change the name of the ribosome data path according to the experiment name used before.**

**Note:** It is recommended to duplicate and alter `template.json` in line with the ribosome experiments applied. Additionally, it is recommended to alter the names of the trained models (defined under `--name` in `train.sh`) in accordance with the ribosome data sets applied.

`template.json`:
```
{
  "h5_path":"../data/GRCh38_v107.h5",
  "exp_path":"transcript",
  "y_path":"tis",
  "chrom_path":"contig",
  "id_path":"id",
  "seq":false,                                  # No sequence data is used as input
  "ribo":{
    experiment_1/5: {}                          <-- Make sure to change this
      }
}
```

To train the models in sequence, run:
```
bash train.sh
```

## ✔️ Roadmap

- [x] Process transcriptome features
- [x] Process ribosome profiling data
- [x] Set-up data format for model training/prediction
- [ ] Post-processing features
- [ ] Model output aggregation and processing 
- [ ] Wrap it all together: 
    - [ ] Single pip package
    - [ ] Simplified/intuitive utility
    - [ ] End-to-end pipeline
- [ ] Optional: support for GUI (streamlit)
