avTAD
=======

DEPRECATION WARNING: This tool was developed before the release of other software for the same purpose (see [coolpuppy](https://github.com/Phlya/coolpuppy)). avTAD uses the deprecated library [mirnylib](https://github.com/mirnylab/mirnylib-legacy), which is no longer maintained, but still can be used with this package. 

avTAD is a tool for construction of average TAD plots. Basic pipeline consists of three steps: *snip*, *rescale* and *plot*. 

### Tool structure

avTAD has four subcommands:
- **snip** - reads TSV (BED-like) file with TADs and input Hi-C map in COOL or HICLIB format; produces TSV table with 
TADs contacts enrichment (log2 of observed over expected) and PICKLE file with individual snips.
- **rescale** - reads PICKLE file wih snips, rescales them to the same size and calculates averaging function over them.
Result is written as TSV file.
- **evaluate** - pipeline extension for comparison of average TADs. Evaluation of simple operations like difference between average TADs of the same size prodiced by *rescale*. 
- **plot** - reads TSV file with average TAD matrix and plots a heatmap. TSV matrix file can be either *rescale* or *evaluate* output. 

### Installation and requirements

Requirements:
- Linux-based OS
- Python >= 3.5

Installation with pip:

```bash
git clone https://github.com/agalitsyna/avTAD.git
cd avTAD
pip install -r requirements.txt
pip install https://bitbucket.org/mirnylab/mirnylib/get/tip.tar.gz
pip install -e .
```

Required packages as listed in requirements.txt (automatically installed with pip install requirements.txt): 
- hdf5>=1.10.4 & h5py>=2.9.0
- cooler>=0.8.5
- numpy>=1.16.4
- pandas>=0.24.2
- matplotlib>=3.1.0
- seaborn>=0.9.0
- click>=7.0
- pytest>=4.6.2
- mirnylib from https://bitbucket.org/mirnylab/mirnylib/get/tip.tar.gz
- biopython (for mirnylib)
- joblib>=0.6.3 (for mirnylib)
- Cython>=0.16 (for mirnylib)


### Example runs
```bash
avTAD snip data/OSC_TADS.bed data/OSC_dm3.cool tmp_results --format cool --diagonals-to-remove 2 --niter 2
avTAD rescale OSC.TADsnips.pickle OSC.TADmetadata.tsv OSC --rescaled-size 200 
avTAD plot OSC OSC --vmax 0.1
```

Example split by chromosomes:
```bash
avTAD rescale OSC.TADsnips.pickle OSC.TADmetadata.tsv OSC --split-by ch --rescaled-size 200 
avTAD plot OSC OSC --vmax 0.1
```

Example input from stdin, grep one chromosome:
```bash
avTAD rescale OSC.TADsnips.pickle <(grep chrX OSC.TADmetadata.tsv) OSC_chrX  --rescaled-size 200 
avTAD plot OSC_chrX OSC_chrX --vmax 0.1
```

Example query:
```bash
avTAD rescale OSC.TADsnips.pickle OSC.TADmetadata.tsv OSC_TADsize10-30  --rescaled-size 200 --query "TAD_size>10 and TAD_size<30"
avTAD plot OSC_TADsize20-30 OSC_TADsize10-30 --vmax 0.1
```

Example comparison with shuffle:
```bash
avTAD snip data/OSC_TADS.bed data/OSC_dm3.cool OSC --format cool --diagonals-to-remove 2 --niter 1
avTAD rescale OSC.TADsnips.pickle OSC.TADmetadata.tsv OSC --rescaled-size 200 
avTAD rescale OSC.TADsnips_shuf0.pickle OSC.TADmetadata.tsv OSC_shuf0  --rescaled-size 200
avTAD evaluate OSC OSC_shuf0 OSC_enrichment "a-b"
avTAD plot OSC OSC --vmax 0.1
avTAD plot OSC_shuf0 OSC_shuf0 --vmax 0.1
avTAD plot OSC_enrichment OSC_enrichment --autoscale
```

Example comparison between experiments:
```bash
avTAD snip data/OSC_TADS.bed data/OSC_dm3.cool OSC --format cool --diagonals-to-remove 2 --niter 1
avTAD snip data/BG3_TADS.bed data/BG3_dm3.cool BG3 --format cool --diagonals-to-remove 2 --niter 1
avTAD rescale OSC.TADsnips.pickle OSC.TADmetadata.tsv OSC --rescaled-size 200 
avTAD rescale BG3.TADsnips.pickle BG3.TADmetadata.tsv BG3 --rescaled-size 200 
avTAD evaluate OSC BG3 OSC-BG3 "a-b"
avTAD plot OSC OSC --vmax 0.1
avTAD plot BG3 BG3 --vmax 0.1
avTAD plot OSC-BG3 OSC-BG3 --autoscale
```
