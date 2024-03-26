<div align="center">
<h1>🧮 RiboTIE </h1>

*Driving coding sequence discovery since 2023*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8059764.svg)](https://doi.org/10.5281/zenodo.8059764)
[![GitHub license](https://img.shields.io/github/license/jdcla/RIBO_former)](https://github.com/jdcla/RIBO_former/blob/main/LICENSE.md)
[![GitHub issues](https://img.shields.io/github/issues/jdcla/RIBO_former)](https://github.com/jdcla/RIBO_former/issues)
[![GitHub stars](https://img.shields.io/github/stars/jdcla/RIBO_former)](https://github.com/jdcla/RIBO_former/stargazers)


</div>

## 📋 About
**Note that RiboTIE was formerly known as RIBO-former.** Changes in naming are ongoing.

[RiboTIE](https://doi.org/10.1101/2023.06.20.545724) is created to annotate translation initiation sites on transcripts using ribosome profiling data. This repository contains the instructions to run RiboTIE on custom data.

The data, model parameters, and benchmark data from **the article** are featured in a [separate repository](https://github.com/jdcla/RIBO_former_paper).

When interested in more advanced features, such as using a custom transformer architecture, we refer the user manual of the [transcript-transformer package](https://github.com/jdcla/transcript_transformer),  created in support of RiboTIE. 

Make sure to check out [TIS Transformer](https://github.com/jdcla/TIS_transformer) as well, a similar tool for the delineation of novel coding sequences using transcript sequence data rather than ribosome profiling data.



## 📖 User guide

### System Requirements
Following are the instructions on how to set up RiboTIE and pre-process data. RiboTIE is a deep learning framework that requires use of a single GPU with CUDA support . All software runs on packages installed through Python

### Installation

`PyTorch` is used as the deep learning library. Follow the instructions [here](https://pytorch.org/get-started/locally/) to install `PyTorch` first. GPU support is necessary.

After installing `PyTorch`, install the `transcript_transformer` packages through PyPI by running

```bash
pip install transcript_transformer
```

Installation should take no longer than a couple of minutes.

### Usage

Before running the tool, a dictionary file (yaml) needs to exist that points towards all the input data used. In addition, this file specifies how data is used to train and evaluate riboformer models. Inspect `template.yml` to evaluate all available options. see `test/` for example inputs. Required are:

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

A pre-trained model is used as this improves performances while drastically reducing computational resources required for the fine-tuning as compared to training models from scratch.

To run RiboTIE:
```bash
riboformer template.yml
```

For more information about various other options, try:
```bash
riboformer -h
```

### Parsing data
Parsing data can be achieved without doing fine-tuning and prediction by running:

```bash
riboformer template.yml --data
```

Once completed, the tool will automatically skip to the fine-tuning and prediction steps when re-running the script (i.e., `riboformer yaml_file.yml`), as the data is detected within the `h5` database.

### Results

RiboTIE evaluates and returns all positions on the transcriptome (saved in `*.npy` files). In addition, RiboTIE collects metadata for the top results within a result table (both `*.csv` and `*.gtf`) for further evaluation. Different flags are available to filter down the results in the output table:

- `--prob_cutoff` (default=0.15) : The model output threshold with which to determine the positive set.
- `--start_codons` (default="\*TG") : Regular expression denoting valid start codons. If all start are viable, use "*".
- `--min_ORF_len` (default=15) : Minimum nucleotide length of predicted translated ORF.
- `--include_invalid_TTS` (default=False) : Include ORFs with no valid stop codon.

The default parameters are our recommendations. To adjust the outputs after having run RiboTIE, make sure to re-run the code with the `--results` flag to prevent the software from re-processing the samples from scratch. For more steps, we include a plethora of metadata in the output tables that can be used to filter against (e.g., `ORF_type`, `tr_support_lvl`, `tr_biotype`, ...). 

In addition to some basic filtering of ORFs, sites with near-miss predictions are corrected ([explanation](https://github.com/jdcla/RIBO_former/blob/main/README.md#near-miss-identifier)). 

Example: create result tables without applying near-miss correction:

```bash
riboformer yaml_file.yml --results --no-correction
```

## pre-trained models

Currently, only a single set of pre-trained models is available: `50perc_06_23.ckpt`. This model is selected by default.
The models are pre-trained on non-overlapping parts of the transcriptome, using the data featured by SRR592960, SRR1562539, SRR1573939, SRR1610244, SRR1976443, SRR2536856, SRR2873532, SRR3575904.

|        |  Train                  | Validation | Applied on                    |
|--------|------------------------|------------|--------------------------------|
| Fold 1 |  3, 5, 7, 11, 13, 15, 19, 21, X | 1, 9, 17     | 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, Y, ... |
| Fold 2 |  2, 6, 8, 10, 14, 16, 18, 22, Y | 4, 12, 20    | 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, X  |
Listed identifiers refer to chromosomes.

Adding the ability to create custom pre-trained models through the `riboformer` script is planned for the future. This can already be achieved when running your full scripts through functionalities within the `transcript-transformer` package.


## How does RiboTIE work?

See [the manuscript](https://www.biorxiv.org/content/10.1101/2023.06.20.545724v1) for a detailed description.
Essentially, RiboTIE works by detecting translation initiation sites using only ribosome profiling data. The tool parses information on how reads are aligned along the transcript. Specifically, for each position, a vector containing the number of reads for each read length at that position is parsed.  **No sequence information is processed**. RiboTIE similarly returns predictions for each position on each transcript.

<div align="center">
<img src="https://github.com/jdcla/RIBO_former/raw/main/ribo_intro.png" width="800">
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

## How can RiboTIE improve?

This list is non-exhaustive, but rather lists low-hanging fruit. 

### Calibration

In line with good machine learning practice, models are not used to obtain predictions on data it is trained on. RiboTIE therefore trains/fine-tunes multiple models on non-overlapping folds of the transcriptome. Predictions over the full transcriptome are gathered by simply merging the outputs of both models. As post-processing and filtering of sites of interest is done on a rank-based merit, this technique is not optimal. In other words, the output distributions are not necessarily aligned, where an output of 0.6 for one model is *as significant* as a 0.6 for the other model (for two-fold approaches).

**Objective**: Apply calibration steps that seeks to improve the creation of a merged ranking when combining multiple sets of predictions from different folds of the data.

**Note:** evaluating PR/ROC AUC of the independent sets and the combined set showed only a very slight decline in performance (~1%). At this point, no time has been invested to address this issue.

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
