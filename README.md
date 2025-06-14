<div align="center">
<h1>🧮 RiboTIE </h1>

*Driving coding sequence discovery since 2023*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8059764.svg)](https://doi.org/10.5281/zenodo.8059764)
[![GitHub license](https://img.shields.io/github/license/TRISTAN-ORF/RiboTIE)](https://github.com/TRISTAN-ORF/RiboTIE/blob/main/LICENSE.md)
[![GitHub issues](https://img.shields.io/github/issues/TRISTAN-ORF/RiboTIE)](https://github.com/TRISTAN-ORF/RiboTIE/issues)
[![GitHub stars](https://img.shields.io/github/stars/TRISTAN-ORF/RiboTIE)](https://github.com/TRISTAN-ORF/RiboTIE/stargazers)


</div>

## 📋 About

[RiboTIE](http://biorxiv.org/cgi/content/full/2024.03.21.586110v1) is created to find translated ORFs on transcripts using ribosome profiling data. RiboTIE achieves state-of-the-art results independent of the lab protocol (e.g., translation blockers) applied to create the ribosome profiling database. This repository contains the instructions to run RiboTIE on custom data.

The data, model parameters, and benchmark data from **the article** are featured in a [separate repository](https://github.com/TRISTAN-ORF/RiboTIE_article).

When interested in more advanced features, such as using a custom transformer architecture, we refer the user manual of the [transcript-transformer package](https://github.com/TRISTAN-ORF/transcript_transformer),  created in support of RiboTIE. 

Make sure to check out [TIS Transformer](https://github.com/TRISTAN-ORF/TIS_transformer) as well, a similar tool for the delineation of novel coding sequences using transcript sequence data rather than ribosome profiling data. TIS Transformer predicts novel CDSs based on transcript sequence data rather than ribosome profiling data, and can be used in conjuction with RiboTIE.



## 📖 User guide

Following are the instructions on how to set up RiboTIE and pre-process data.

### Installation

`PyTorch` is used as the deep learning library. Follow the instructions [here](https://pytorch.org/get-started/locally/) to install `PyTorch` first. GPU support is necessary.

After installing `PyTorch`, run

```bash
pip install transcript_transformer
```

### Usage

**All commands listed here can be executed within this repository.**

Dictionary files (YAML/JSON) are the recommended approach to pass arguments to the tool. It is possible to list multiple configuration files. Inspect `template.yml` and `ribotie -h` to evaluate all available options. see `test/` for example inputs. Required are:

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
5. Filter out [CDS variant calls](https://github.com/TRISTAN-ORF/RiboTIE/blob/main/README.md#variant-cds-filtering)

RiboTIE finetunes a pre-trained model on individual samples as this improves performances while drastically reducing computational resources required during the fine-tuning step as compared to training models from scratch. By default, RiboTIE incorporates a pre-trained model for human data. Alternatively, [models can be pre-trained on custom data](https://github.com/TRISTAN-ORF/RiboTIE/blob/main/README.md#pre-training-models), e.g., to apply RiboTIE on other species. 

To run RiboTIE on some test data, install `transcript_transformer`, clone this repository and run:
```bash
ribotie template.yml
```

For more information about various other options, try:
```bash
ribotie -h
```

### Parsing data
Parsing the data within the `hdf5` database can be achieved without doing fine-tuning and prediction. This can be specifically of interest when running RiboTIE on cloud infrastructure.

```bash
ribotie template.yml --data
```

Once completed, the tool will automatically skip to the fine-tuning and prediction steps when re-running the script (i.e., `ribotie template.yml`) as the data is detected within the `hdf5` database.

### Parsing the predictions 

RiboTIE evaluates and returns all positions on the transcriptome (saved in `*.npy` files). It is not feasible or of interest to provide information on the millions of predictions available.
Therefore, RiboTIE only collects metadata for predictions meeting the [listed criteria](https://github.com/TRISTAN-ORF/RiboTIE/blob/main/README.md#filtering). For these, the tool will generate a result table (`*.csv`) and a minimally formatted `.gtf` file that can be combined with tools such as `gffread` to extract sequences. 

**NOTE: the output `.csv`/`.gff` files are generated from the full set of predictions on the transcriptome (within the `.npy` files). The RiboTIE model should not be fine-tuned again when creating new result tables (Use the `--results` flag!)** 

It is possible to alter the conditions to select called ORFs when generating the metadata table.

For example: get more metadata on the ORFs where all start codons are included:

```bash
ribotie yaml_file.yml --results --start_codons ".*"
```


### Evaluating results

RiboTIE now generates several report files that are automatically detected by MultiQC and displayed when running it in a parent directory.

### Pre-trained models


Currently, a single set of pre-trained models is available which can be used for human ribosome profiling data. This model is the selected one by default.
The models are pre-trained on non-overlapping parts of the transcriptome, using the data featured by SRR592960, SRR1562539, SRR1573939, SRR1610244, SRR1976443, SRR2536856, SRR2873532, SRR3575904.

|        |  Train                  | Validation | Applied on                    |
|--------|------------------------|------------|--------------------------------|
| Fold 1 |  3, 5, 7, 11, 13, 15, 19, 21, X | 1, 9, 17     | 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, Y, ... |
| Fold 2 |  2, 6, 8, 10, 14, 16, 18, 22, Y | 4, 12, 20    | 1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, X  |
Listed identifiers refer to chromosomes.

### Pre-training models 

It is possible to pre-train custom models. This is highly recommended when pre-trained models are not available (e.g., new organisms). If your study includes multiple samples, it is possible pre-train models on all samples before the fine-tuning step.  Note that two models are optimized to ensure no unforeseen biases are introduced by training or selecting (validation set) models on identical genomic regions they are providing predictions for (also dubbed test set). 

```bash
ribotie template.yml --pretrain
```

The following steps are performed:
1. RiboTIE divides the transcripts by seqname into training, validation and test sets using two folds. The test sets cover the full transcriptome (50% in each fold). The test set is not used during training or selection of the final models.
2. RiboTIE combines all listed samples (`--samples`) into a single training, validation and test set for each fold/model.
3. After optimization, RiboTIE saves the weights of the models (`.ckpt`) as `{out_prefix}pretrain_f1.ckpt` and `{out_prefix}pretrain_f2.ckpt`. It also creates a `yaml` output file (under `{out_prefix}pretrain.yml`), which contains all ribotie arguments required to apply these newly pretrained models, as well as the genomic regions allocated for generating the training, validation and test sets.

To run RiboTIE using custom pre-trained models:

```bash
ribotie template.yml out_test/pretrain.yml
```

Note that RiboTIE supports using multiple separate configuration files (the order is not important). By default, RiboTIE will determine the directory path of the configuration file listing the pre-trained model to search for the model dictionary (`.ckpt`). Both files can be moved together to a new directory if desired. 

It is furthermore possible to determine the allocations for the training, validation and test sets manually.

```bash
ribotie template.yml test/folds.yml --pretrain
```


### ⚡️ Parallelization

The use of a single `hdf5` database has limitations towards upscaling and parallelization of RiboTIE. While samples are processed independently, having multiple processes write and read to and from a single `hdf5` can result in I/O errors.

The argument `--parallel` creates independent `hdf5` datasets for each ribo-seq sample  (`{h5_path%.h5}_{ribo_id}.h5`). This allows individual samples being processed in parallel by RiboTIE. RiboTIE still pulls data from the `hdf5` database containing the genomic features when processing the result table, so make sure these databases always exist under the same folder as defined by the `--h5_path` argument.

For example, a snakemake pipeline parsing `ribo-seq` samples in parallel using RiboTIE could look like this:

```python
# from template.yml, can also be part of a snakemake configuration file
samples = {
    "sample_1": "test/SRR000001.sam",
    "sample_2": "test/SRR000001.sam",
    "sample_3": "test/SRR000001.sam"
}

rule all:
    input:
        expand("test_out/{sample}.npy", sample=samples.keys())

# `tis_transformer --data` is called first as running multiple `ribotie` commands in parallel
# would result in parallel processes parsing the genome assembly features, which is only required
# once (and stored under `{h5_path}` and within the backup location).
rule ribotie_parse_genomic_features:
    input:
        "template.yml"
    output:
        "GRCh38v110_snippet.h5"
    shell:
        """
        tis_transformer {input} --data
        """

# For HPC servers, it can further be advantageous to create an extra rule for parsing the riboseq
# data separately using a partition that does not require GPU resources
# e.g. `ribotie {input.config} --data --samples {wildcards.sample} --parallel`
rule ribotie_parse_and_predict_riboseq_samples:
    input:
        config="template.yml",
        base="GRCh38v110_snippet.h5",
        mapped=lambda wildcards: samples[wildcards.sample]
    output:
        "test_out/{sample}.npy"
    shell:
        """
        ribotie {input.config} --samples {wildcards.sample} --parallel
        """

```

> Note that if a variable is defined both within yaml configuration file and a bash argument (e.g. `--samples`), the latter (i.e., bash arguments) will overwrite the former.  

### Filtering

#### Custom filters

Custom filtered can be toggled or altered by the user. By default, included predictions are filtered on:

- `--prob_cutoff`: the model prediction score (`ribotie_score`) is larger than 0.125
- `--start_codons`: start codons (`start_codon`) are near-cognate (*TG)
- `--include_invalid_TTS`: a valid translation termination site (`TTS_on_transcript`) is present on the transcript.

#### Variant CDS filtering

<div align="center">
<img src="https://github.com/TRISTAN-ORF/RiboTIE/raw/main/filtering.png" width="800">
</div>

Results excluding variant CDS filtering are generated as `*.unfiltered.gtf` and `*.unfiltered.csv`. Both the custom filters and near-miss correction can be toggled through the command line tool.

#### Near-miss identifier

RiboTIE, unlike previous tools processing ribosome profiling data, does not create ORF libraries or has access to start codon information when making predictions. Essentially, it only parses ribosome profiling information along the transcript.

It is observed that, for transcripts featuring fewer mapped reads around the translation initiation site, RiboTIE is more prone to miss translation initiation sites by several bases. To address this issue, a neighborhood searching step is performed when creating the result table that corrects **non-ATG** predictions to **in-frame ATG positions**  if **present within a 9 codons distance**. Performed corrections are listed as `correction` in the result table. This feature can be disabled using the `--no-correction` flag. 

## 🤨 How does RiboTIE work?

See [the manuscript](https://www.biorxiv.org/content/10.1101/2023.06.20.545724v1) for a detailed description.
RiboTIE detects translation initiation sites using only ribosome profiling data. The tool parses information on how reads are aligned along the transcript. Specifically, for each position, a vector containing the number of reads for each read length at that position is parsed.  **No sequence information is processed**. RiboTIE similarly returns predictions for each position on each transcript. From these predictions, ORFs are derived by a greedy search algorithm that finds the first in-frame stop codon. 

<div align="center">
<img src="https://github.com/TRISTAN-ORF/RiboTIE/raw/main/ribo_intro.png" width="800">
</div>

**Note:** striked-through text refers to steps typically performed by existing methods but ommited by RiboTIE.

Fine-tuning is important as ribosome profiling data has shown to be highly variable between experiments, with biological variability (tissues, sample age, ...) and technical variability (lab protocol, machine calibration, ...) playing a role.
RiboTIE is trained and fine-tuned using a set of canonical coding sequences. This approach does not prevent the trained model to find non-canonical ORFs.
After fine-tuning, the model provides predictions for each  position of each transcript. 
No additional post-processing steps are performed.
From these predicted translation initiation sites, the resulting translated ORFs are obtained by searching for the first in-frame stop codon.
Filtering can be applied based on the characteristics of the translated ORFs (e.g. start codon, minimum length).

This technique was shown to substantially outperform previous methods. We hypothesize this gain to be achieved through various factors:
- fine-tuning on each data set, the model learns custom rules present for each data set
- inclusion of read length information
- elegant approach with very few custom hardcoded rules for data (pre-)processing or selection.
- use of a state-of-the-art machine learning tools (transformer networks), which are perfectly suited for the data type (a variable number of input vectors). 


## (optional) Post-processing steps that correct RiboTIE predictions

## ⁉️ FAQ

**Will training and selecting models on canonical coding sequences limit the ability of RiboTIE to detect non-canonical ORFs?**

RiboTIE has been trained only on ribosome sequencing data to prevent the model from relying on ORF metadata or codon information, from which the length of ORFs is easily derived. 
As such, RiboTIE would only be biased against detecting TISs of shorter ORFs if the distribution of ribosome protected fragments is distinctly different from that of canonincal coding sequences, for which there exists no proof today. 
This statement is reflected by our benchmark of RiboTIE against other existing tools, where the discovery rate of short canonical CDSs (<300nt) by RiboTIE was at least 400% higher than any other tool [main manuscript; Fig 1C](http://biorxiv.org/cgi/content/full/2024.03.21.586110v1). 

## 🖊️ Citation  

```bibtex
@article{clauwaert2025deep,
  title={Deep learning to decode sites of RNA translation in normal and cancerous tissues},
  author={Clauwaert, Jim and McVey, Zahra and Gupta, Ramneek and Yannuzzi, Ian and Basrur, Venkatesha and Nesvizhskii, Alexey I and Menschaert, Gerben and Prensner, John R},
  journal={Nature Communications},
  volume={16},
  number={1},
  pages={1275},
  year={2025},
  publisher={Nature Publishing Group UK London}
}
```

## ✔️ Roadmap

- [x] Process transcriptome features
- [x] Process ribosome profiling data
- [x] Set-up data format for model training/prediction
- [x] Post-processing features
    - [x] Result table (top predictions)
    - [x] Assess near-miss predictions
- [ ] Usability
    - [x] User-defined filtering
    - [ ] User-defined output data
        - [x] csv file (containing called ORF/transcript/gene metadata)
        - [x] gtf file (containing called ORF genomic features)
        - [ ] fasta file (containing called ORFs sequences)
    - [x] Support pre-training models on custom sets of data
    - [x] Parallelization for use with NextFlow/Snakemake
- [ ] Optimizations
    - [ ] SAM/BAM parser that does not load file completely in memory
- [x] Wrap it: 
    - [x] Single pip package
    - [x] Simplify README
    - [x] Create a custom call function
    - [x] Add pre-trained models to package, pre-define required input arguments
    - [x] End-to-end pipeline
- [ ] Optional: support for GUI (streamlit)
