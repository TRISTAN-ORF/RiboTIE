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

When interested in more advanced features, such as using a custom transformer architecture, we refer the user manual of the [transcript-transformer package](https://github.com/jdcla/transcript_transformer),  created in support of RIBO-former. 

Make sure to check out [TIS Transformer](https://github.com/jdcla/TIS_transformer) as well, a similar tool for the delineation of novel coding sequences using transcript sequence data rather than ribosome profiling data.



## 📖 User guide

Following are the instructions on how to set up RIBO-former and pre-process data.

### Installation

`PyTorch` is used as the deep learning library. Follow the instructions [here](https://pytorch.org/get-started/locally/) to install `PyTorch` first. GPU support is necessary.

After installing `PyTorch`, run

```bash
pip install transcript_transformer
```

### Usage

Before running the tool, a dictionary file (YAML) needs to exist that points towards all the input data used. In addition, this file specifies how data is used to train and evaluate riboformer models. Inspect `template.yml` to evaluate all available options. Required are:

- reference assembly files (`*.gtf`, `*.fa`)
- ribosome profiling reads (`*.sam`) **mapped to the transcriptome**

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

When running RIBO-former, the following steps are performed:

1. Parse all data to a HDF5 database (`h5_path`)
2. Fine-tune pre-trained models on non-overlapping folds of the data. This allows the model to learn data-set specific correlations that are relevant.
3. Get model predictions for all positions of the transcriptome
4. Collect metadata for the top ranking predictions

A pre-trained model is used as this improves performances while drastically reducing computational resources required for the fine-tuning as compared to training models from scratch.

To run RIBO-former:
```bash
riboformer yaml_file.yml
```

For more information about specific options, try:
```bash
riboformer -h
```

**Note**: Currently, parsing data from `.sam` files to the `h5` database can require high amounts of RAM. For multiple data sets, this process might run for several hours. It is possible to pre-process the data without doing fine-tuning and prediction by running:

```bash
riboformer yaml_file.yml --data-process
```
Afterwards, running `riboformer yaml_file.yml`, the module will automatically detect data already present in the `h5` database, skipping to the fine-tuning step.

## pre-trained models

Currently, only a single set of pre-trained models is available: `50perc_06_23.ckpt`. This model is selected by default.
The models are pre-trained on non-overlapping parts of the transcriptome, using the data featured by SRR592960, SRR1562539, SRR1573939, SRR1610244, SRR1976443, SRR2536856, SRR2873532, SRR3575904.

|        |  Train                  | Validation | Applied on                    |
|--------|------------------------|------------|--------------------------------|
| Fold 1 |  3, 5, 7, 11, 13, 15, 19, 21, X | 1, 9, 17     | 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, Y, ... |
| Fold 2 |  2, 6, 8, 10, 14, 16, 18, 22, Y | 4, 12, 20    | 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, X  |
Listed identifiers refer to chromosomes.

Adding the ability to create custom pre-trained models through the `riboformer` script is planned for the future. This can already be achieved when running your full scripts through functionalities within the `transcript-transformer` package.


## How does RIBO-former work?

See [the manuscript](https://www.biorxiv.org/content/10.1101/2023.06.20.545724v1) for a detailed description.
Essentially, RIBO-former works by detecting translation initiation sites using only ribosome profiling data using transformer networks. The tool parses information on how reads are aligned along the transcript. Specifically, for each position, a vector containing the number of reads for each read length at that position is parsed.  **No sequence information is processed**. Ribo-former similarly returns predictions for each position on each transcript.

<div align="center">
<img src="https://github.com/jdcla/RIBO_former/raw/main/ribo_intro.png" width="800">
</div>

**Note:** striked-through text refers to steps typically performed by existing methods but ommited by RIBO-former.

Fine-tuning is important as ribosome profiling data has shown to be highly variable between experiments, with biological variability (tissues, sample age, ...) and technical variability (lab protocol, machine calibration, ...) playing a role.
RIBO-former is trained and fine-tuned using a set of canonical coding sequences. This approach does not prevent the trained model to find non-canonical ORFs.
The script simply returns the top ranking predictions of the full set of predictions evaluated on each position of each transcript. 
No additional post-processing steps are performed.
From these predicted translation initiation sites, the resulting translated ORFs are obtained by searching for the first in-frame stop codon.
No filtering is applied based on the characteristics of the translated ORFs (e.g. start codon, minimum length).

This technique was shown to substantially outperform previous methods. We hypothesize this gain to be achieved through various factors:
- fine-tuning on each data set, the model learns custom rules present for each data set
- inclusion of read length information
- elegant approach with very few custom rules for data (pre-)processing or selection.

## How can RIBO-former improve?

This list is non-exhaustive, but rather lists low-hanging fruit. 

### Calibration

In line with good machine learning practice, models are not used to obtain predictions on data it is trained on. RIBO-former therefore trains/fine-tunes multiple models on non-overlapping folds of the transcriptome. Predictions over the full transcriptome are gathered by simply merging the outputs of both models. As post-processing and filtering of sites of interest is done on a rank-based merit, this technique is not optimal. In other words, the output distributions are not necessarily aligned where an output of 0.6 for one model is *as significant* as a 0.6 for the other model (for two-fold approaches).

**Objective**: Apply calibration steps that seeks to improve the ranking of multiple sets of predictions from different folds of the data.

Note: evaluating PR/ROC AUC of the independent sets and the combined set showed only a slight decline in performance (~2%). As such, no time has been invested to implement improvements at this point.

### Near-miss identifier

RIBO-former, unlike previous tools processing ribosome profiling data, does not create ORF libraries or has access to start codon information. Essentially, it only parses ribosome profiling information along the transcript.

It is observed that, for transcripts featuring a lower number of mapped reads (low coverage and read depth), RIBO-former can miss the exact location of well-known translation initiation sites by several bases.

**Objective**: Implement a neighborhood searching step that evaluates near-miss predictions when processing a set of top-ranking predictions. This information can be included as part of the metadata.


## ✔️ Roadmap

- [x] Process transcriptome features
- [x] Process ribosome profiling data
- [x] Set-up data format for model training/prediction
- [ ] Post-processing features
    - [X] Result table (top predictions)
    - [ ] Calibrate predictions from different folds/models
    - [ ] Assess near-miss predictions
- [ ] Usability
    - [ ] User-defined filtering
    - [ ] User-defined output formatting
    - [ ] Pre-training models on custom sets of data
- [x] Wrap it: 
    - [x] Single pip package
    - [x] Simplify README
    - [x] Create a custom call function
    - [x] Add pre-trained models to package, pre-define required input arguments
    - [x] End-to-end pipeline
- [ ] Optional: support for GUI (streamlit)
