<div align="center">
<h1>🧮 RiboTIE </h1>

*Driving coding sequence discovery since 2023*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8059764.svg)](https://doi.org/10.5281/zenodo.8059764)
[![GitHub license](https://img.shields.io/github/license/TRISTAN-ORF/RiboTIE)](https://github.com/TRISTAN-ORF/RiboTIE/blob/main/LICENSE.md)
[![GitHub issues](https://img.shields.io/github/issues/TRISTAN-ORF/RiboTIE)](https://github.com/TRISTAN-ORF/RiboTIE/issues)
[![GitHub stars](https://img.shields.io/github/stars/TRISTAN-ORF/RiboTIE)](https://github.com/TRISTAN-ORF/RiboTIE/stargazers)


</div>

## 📋 About
**Note that RiboTIE was formerly known as RIBO-former.** Changes in naming are ongoing.

[RiboTIE](https://doi.org/10.1101/2023.06.20.545724) is created to annotate translation initiation sites on transcripts using ribosome profiling data. This repository contains the instructions to run RiboTIE on custom data.

The data, model parameters, and benchmark data from **the article** are featured in a [separate repository](https://github.com/TRISTAN-ORF/RiboTIE_article).

When interested in more advanced features, such as using a custom transformer architecture, we refer the user manual of the [transcript-transformer package](https://github.com/TRISTAN-ORF/transcript_transformer),  created in support of RiboTIE. 

Make sure to check out [TIS Transformer](https://github.com/TRISTAN-ORF/TIS_transformer) as well, a similar tool for the delineation of novel coding sequences using transcript sequence data rather than ribosome profiling data.



## 📖 User guide

Following are the instructions on how to set up RiboTIE and pre-process data.

### Installation

`PyTorch` is used as the deep learning library. Follow the instructions [here](https://pytorch.org/get-started/locally/) to install `PyTorch` first. GPU support is necessary.

After installing `PyTorch`, run

```bash
pip install transcript_transformer
```

### Usage

> Note that all commands listed here can be executed within the directory after cloning the repository

Dictionary files (YAML/JSON) are the recommended approach to pass arguments to the tool. It is possible to list multiple configuration files. Inspect `template.yml` to evaluate all available options. see `test/` for example inputs. Required are:

- a **genome-level** reference and assembly file (`*.gtf`, `*.fa`)
- ribosome profiling reads (`*.sam`/`*.bam`) **mapped to the transcriptome**

```yaml
gtf_path : path/to/gtf_file.gtf
fa_path : path/to/fa_file.fa
########################################################
## add entries to ribosome profiling data.
## format: 'id : ribosome profiling paths'
########################################################
ribo_paths :
  SRR000001 : path/to/mapped/riboseq.sam
  SRR000002 : path/to/mapped/riboseq.sam
  SRR000003 : path/to/mapped/riboseq.sam
########################################################
## database path (parsed data output)
########################################################
h5_path : my_experiment.h5
```

When running RiboTIE, the following steps are performed:

1. Parse all data to a HDF5 database (`h5_path`)

subsequently, for every data set in `ribo_paths`:

2. Fine-tune pre-trained models on non-overlapping folds of the data.
3. Get model predictions for all positions of the transcriptome
4. Collect metadata for the top ranking predictions

RiboTIE finetunes a pre-trained model on individual samples as this was shown to improve performances while drastically reducing computational resources required during the fine-tuning step as compared to training models from scratch. By default, RiboTIE incorporates a pre-trained model for human data. See "Pre-training a custom model" for instructions on how to pre-train a custom model. e.g., to apply RiboTIE on other species. 

To run RiboTIE on some test data, clone this repository and run:
```bash
ribotie template.yml
```

For more information about various other options, try:
```bash
ribotie -h
```

### Parsing data
Parsing data can be achieved without doing fine-tuning and prediction by running:

```bash
ribotie template.yml --data
```

Once completed, the tool will automatically skip to the fine-tuning and prediction steps when re-running the script (i.e., `ribotie yaml_file.yml`), as the data is detected within the `hdf5` database.

### Parsing the predictions 

RiboTIE evaluates and returns all positions on the transcriptome (saved in `*.npy` files). It is not feasible or of interest to provide information on the millions of predictions.
Therefore, RiboTIE collects metadata for predictions meeting several criteria. For these, the tool will generate a result table (`*.csv`) and a minimally formatted `.gtf` file that can be combined with tools such as `gffread` to extract sequences. 

By default, included predictions are filtered on:

- the model output (`output`) is larger than 0.15
- start codons (`start_codon`) are near-cognate (*TG)
- a valid translation termination site (`TTS_on_transcript`) is present on transcript

Additionally, sites with near-miss predictions are corrected ([explanation](https://github.com/TRISTAN-ORF/RiboTIE/blob/main/README.md#near-miss-identifier)).

**NOTE: the output `.csv`/`.gff` files are generated from the full set of predictions on the transcriptome (within the `.npy` files). The RiboTIE model should not be fine-tuned again when creating new result tables (Use the `--results` flag!)** 

It is possible to alter the conditions used to generate the final results.

For example: to include all start codons:

```bash
ribotie yaml_file.yml --results --start_codons ".*"
```

### Evaluating results

RiboTIE know generates several report files that are automatically detected by MultiQC and displayed when running it in a parent directory.

## Pre-trained models


Currently, a single set of pre-trained models is available which can be used for human ribosome profiling data. This model is the selected one by default.
The models are pre-trained on non-overlapping parts of the transcriptome, using the data featured by SRR592960, SRR1562539, SRR1573939, SRR1610244, SRR1976443, SRR2536856, SRR2873532, SRR3575904.

|        |  Train                  | Validation | Applied on                    |
|--------|------------------------|------------|--------------------------------|
| Fold 1 |  3, 5, 7, 11, 13, 15, 19, 21, X | 1, 9, 17     | 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, Y, ... |
| Fold 2 |  2, 6, 8, 10, 14, 16, 18, 22, Y | 4, 12, 20    | 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, X  |
Listed identifiers refer to chromosomes.

Adding the ability to create custom pre-trained models through the `ribotie` script is planned for the future. This can already be achieved when running your full scripts through functionalities within the `transcript-transformer` package.

## Pre-training models 

It is possible to pre-train custom models. This is highly recommended when no pre-trained models are available (e.g., new organisms). If your study includes multiple samples, it is possible pre-train models on all samples before the fine-tuning step.  Note that two models are optimized to ensure no unforeseen biases are introduced by training or selecting (validation set) models on identical genomic regions they are providing predictions for (also dubbed test set). 

```bash
ribotie template.yml --pretrain
```

The following functions are performed:
1. RiboTIE divides the transcripts by seqname into training, validation and test sets using two folds. The test sets cover the full transcriptome (50% in each fold). The test set is not used during training or selection of the final models.
2. RiboTIE combines all listed samples (`--samples`) into a single training, validation and test set for each fold/model.
3. After optimization, RiboTIE saves the weights of the models (`.ckpt`) as `{out_prefix}pretrain_f1.ckpt` and `{out_prefix}pretrain_f2.ckpt`. It also creates a `yaml` output file (under `{out_prefix}pretrain.ckpt`), which contains all ribotie arguments required to apply these newly pretrained models.

To run RiboTIE using custom pre-trained models:

```bash
ribotie template.yml out_test/pretrain.yml
```

Note that RiboTIE supports using multiple separate configuration files (the order is not important). By default, RiboTIE will determine the directory path of the configuration file listing the pre-trained model to search for the model dictionary (`.ckpt`). Both files can be moved together to a new directory if desired. 

## How does RiboTIE work?

See [the manuscript](https://www.biorxiv.org/content/10.1101/2023.06.20.545724v1) for a detailed description.
Essentially, RiboTIE works by detecting translation initiation sites using only ribosome profiling data. The tool parses information on how reads are aligned along the transcript. Specifically, for each position, a vector containing the number of reads for each read length at that position is parsed.  **No sequence information is processed**. RiboTIE similarly returns predictions for each position on each transcript.

<div align="center">
<img src="https://github.com/TRISTAN-ORF/RiboTIE/raw/main/ribo_intro.png" width="800">
</div>

**Note:** striked-through text refers to steps typically performed by existing methods but ommited by RiboTIE.

Fine-tuning is important as ribosome profiling data has shown to be highly variable between experiments, with biological variability (tissues, sample age, ...) and technical variability (lab protocol, machine calibration, ...) playing a role.
RiboTIE is trained and fine-tuned using a set of canonical coding sequences. This approach does not prevent the trained model to find non-canonical ORFs.
The script simply returns the top ranking predictions of the full set of predictions evaluated on each position of each transcript. 
No additional post-processing steps are performed.
From these predicted translation initiation sites, the resulting translated ORFs are obtained by searching for the first in-frame stop codon.
No filtering is applied based on the characteristics of the translated ORFs (e.g. start codon, minimum length).

This technique was shown to substantially outperform previous methods. We hypothesize this gain to be achieved through various factors:
- fine-tuning on each data set, the model learns custom rules present for each data set
- inclusion of read length information
- elegant approach with very few custom hardcoded rules for data (pre-)processing or selection.
- use of a state-of-the-art machine learning tools (transformer networks), which are perfectly suited for the data type (a variable number of input vectors). 

## How RiboTIE predictions are improved

### Near-miss identifier

RiboTIE, unlike previous tools processing ribosome profiling data, does not create ORF libraries or has access to start codon information when making predictions. Essentially, it only parses ribosome profiling information along the transcript.

It is observed that, for transcripts featuring fewer mapped reads around the translation initiation site, RiboTIE is more prone to miss translation initiation sites by several bases. To address this issue, a neighborhood searching step is performed when creating the result table that corrects **non-ATG** predictions to **in-frame ATG positions**  if **present within a 9 codons distance**. Performed corrections are listed as `correction` in the result table. This feature can be disabled by adding the `--no-correction` flag. 

## ✔️ Roadmap

- [x] Process transcriptome features
- [x] Process ribosome profiling data
- [x] Set-up data format for model training/prediction
- [ ] Post-processing features
    - [x] Result table (top predictions)
    - [ ] Calibrate predictions from different folds/models
    - [x] Assess near-miss predictions
- [ ] Usability
    - [ ] Allow pre-trained model on non-human data (detect and split new seqnames evenly in folds)
    - [x] User-defined filtering
    - [ ] User-defined output formatting
    - [ ] Pre-training models on custom sets of data
- [x] Wrap it: 
    - [x] Single pip package
    - [x] Simplify README
    - [x] Create a custom call function
    - [x] Add pre-trained models to package, pre-define required input arguments
    - [x] End-to-end pipeline
- [ ] Optional: support for GUI (streamlit)
