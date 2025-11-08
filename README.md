## Introduction

MGAT (**M**apping-based **G**enome-blocks **A**ncestors **T**racer) is a tool for analysis the ancestors of genome
blocks based on reference genomes.

## Dependencies

### Software

- minimap2
- samtools

### Python Modules

- pathos
- pysam

## Installation

```bash
cd /path/to/install
git clone https://github.com/sc-zhang/MGAT.git
pip install -r MGAT/requirements.txt
chomd +x MGAT/mgat.py

# Optional
echo 'export PATH=/path/to/install/MGAT:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

## Usage

```bash
usage: mgat.py [-h] -r REFERENCE -q QUERY [-w WINDOW] [-s STEP] -o OUTPUT [--minimap_args MINIMAP_ARGS] [--threshold THRESHOLD] [--chain CHAIN] [-t THREADS]

options:
  -h, --help            show this help message and exit
  -r, --reference REFERENCE
                        Directory of reference genomes
  -q, --query QUERY     Directory of query genomes
  -w, --window WINDOW   Window size, default=5000
  -s, --step STEP       Step size, default=5000
  -o, --output OUTPUT   Directory of output files
  --minimap_args MINIMAP_ARGS
                        Arguments for minimap, default="-ax asm10 --eqx"
  --threshold THRESHOLD
                        Threshold for classifying blocks by identity, means best identity should higher than second best identity * threshold, default=1.5
  --chain CHAIN         Block counts for classifying undetermined blocks, the blocks with this counts before current block andthe blocks with this counts after current block would be used, default=5
  -t, --threads THREADS
                        Number of threads, default=10
```

## Result

- 06.results: files in this folder are the final result, which like below:

```text
## Window size: 5000
## Step size: 5000
## Identity threshold: 1.500000
## Chain block counts: 5
#Chrom  Start   Class
Chr10A  1       Ancestor1
Chr10A  5001    Ancestor1
Chr10A  10001   Ancestor2
```