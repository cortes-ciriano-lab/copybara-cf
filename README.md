# ![COPYBARA](/docs/COPYBARA_logo.png)

## Contents
* [Installation](#installation)
  + [Install COPYBARA with Conda](#install-copybara-with-conda)
  + [Install COPYBARA from Source](#install-copybara-from-source)
  + [Check COPYBARA Installation](#check-copybara-installation)
* [Run COPYBARA](#run-copybara)
  + [Panel of Normals Generation](#panel-of-normals-generation)
  + [Focal Amplification Analysis](#focal-amplification-analysis)
  + [Mandatory Arguments](#mandatory-arguments)
  + [Optional Arguments](#optional-arguments)
* [Output Files](#output-files)
* [Troubleshooting](#troubleshooting)
* [License](#license)
* [Contacts](#contacts)

# COPYBARA

COPYBARA is a copy number analysis tool for long-read cell-free DNA (cfDNA) sequencing data. It takes aligned tumour BAM files, performs read counting across genomic bins, normalizes against a panel of normals, matched normal or self, segments the genome using circular binary segmentation (CBS), and fits absolute copy number profiles including tumour fraction, ploidy and percentage of genome abnormality estimations. In addition, COPYBARA focal performs copy number analyses for the detection of focal amplifications and putative extrachromosomal DNA (ecDNA) from cfDNA data. Lastly, COPYBARA also allows _in silico_ size selection during copy number analysis and the generation of custom-made panels of normals (PoN). 

COPYBARA has been tested on Oxford Nanopore (ONT) sequencing reads aligned with minimap2. It requires a Unix-based operating system and has been developed and tested on Linux.

For further information, benchmarking, and citation, please refer to our [preprint](https://github.com/cortes-ciriano-lab/COPYBARA).

## Installation

### Install COPYBARA with Conda

The recommended and easiest way to install COPYBARA is via conda:

```
conda install -c bioconda copybara
```

This will install all dependencies and allow you to use COPYBARA on the command-line.

### Install COPYBARA from Source

**Alternatively**, you can install COPYBARA from source, by first cloning this repository:

```
git clone https://github.com/cmsauer/COPYBARA.git
```

To install from source, COPYBARA requires Python 3.9 with dependencies as listed in the `requirements.txt` file. All dependencies can be installed via conda or pip.

#### Install Dependencies with Conda
To install and manage dependencies with conda, create a new environment and install dependencies (including Python 3.9) with the `environment.yml` file in the top-level of the repository:

```
conda env create --name <env> --file environment.yml
```

#### Install Dependencies with pip
If preferred, you can install and manage dependencies with pip instead using the `requirements.txt` file:

```
pip install -r requirements.txt
```

#### Install COPYBARA
Once you've installed the required dependencies with conda or pip, you can install COPYBARA by navigating to the cloned repo and running:

```
python -m pip install . -vv
```

### Check COPYBARA Installation

You can test that COPYBARA was installed successfully by running `copybara --help`, which should display the following text:

```
usage: copybara [-h] [--version] {pon,focal} ...

COPYBARA - Copy number analysis of long-read cfDNA sequencing data

optional arguments:
  -h, --help     show this help message and exit
  --version      show program's version number and exit

subcommands:
  {pon,focal}
                 COPYBARA sub-commands
    pon          generate panel of normals from bam files for COPYBARA copy number analysis
    focal        detect putative ecDNA/focal amplifications across regions of interest
```

## Run COPYBARA

After installing, COPYBARA can be run on long-read cfDNA data with a minimum set of arguments. The main command to perform copy number analysis with the minimum of required input parameters is:

```
copybara --bam <tumour-bam> --ref <ref-fasta> --outdir <outdir>
```

This will perform read counting, self-normalization, segmentation, and absolute copy number fitting.

To run with a matched normal BAM:

```
copybara --bam <tumour-bam> --normal_bam <normal-bam> --ref <ref-fasta> --outdir <outdir>
```

And, to run with a PoN instead of a matched normal BAM:

```
copybara --bam <tumour-bam> --panel_of_normal <panel-of-normal> --ref <ref-fasta> --outdir <outdir>
```

### Panel of Normals Generation

If you would like to use a panel of normals, for normalisation purposes during copy number analyses, you will need to generate a panel of normals (PoN) from healthy cfDNA samples before running the main copy number command outlined above. The resulting PoN will then provide the baseline read counts for normalization. The PoN can be generated using COPYBARA as follows:

```
copybara pon --pon_list <list-of-normal-bams.txt> --ref <ref-fasta> --outdir <pon-outdir>
```

The `--pon_list` should be a text file with one BAM file path per line.

### Focal Amplification Analysis - COPYBARA focal

COPYBARA also includes a focal analysis mode to detect putative ecDNA or focal amplifications in regions of interest:

```
copybara focal --bam <tumour-bam> --roi <regions-of-interest.bed> --ref <ref-fasta> --outdir <focal-outdir>
```

### COPYBARA Arguments

#### Mandatory Arguments
Argument | Description
-------- | -----------
--bam | Tumour BAM/CRAM file (must have index in .bai/.crai format)
--ref | Full path to reference genome that was used to align the BAM
--outdir | Output directory (can exist but must be empty)

#### Optional Arguments
Argument | Description
-------- | -----------
--normal_bam | Matched normal BAM/CRAM file (must have index in .bai/.crai format)
--panel_of_normal | Path to panel of normal read count file (alternative to --normal_bam)
--sample | Name to prepend to output files (default=tumour BAM filename without extension)
--cn_binsize | Bin window size in kbp (default=500)
--chromosomes | Chromosomes to analyse. To run on all chromosomes leave unspecified (default). To run on a subset of chromosomes only specify the chromosome numbers separated by spaces. For x and y chromosomes, use 23 and 24, respectively. E.g. use "--chromosomes 1 4 23 24" to run chromosomes 1, 4, X and Y
--blacklist | Path to the blacklist file
--blacklist_buffer | Length of region (in bp) flanking the blacklisted region to be excluded (default=0)
--bl_threshold | Percentage overlap between bin and blacklist threshold to tolerate for read counting (default=1)
--bases_threshold | Percentage of known bases per bin required for read counting (default=75)
--smoothing_level | Size of neighbourhood for smoothing (default=10)
--trim | Trimming percentage to be used (default=0.025)
--min_segment_size | Minimum size for a segment to be considered (default=2500000)
--shuffles | Number of permutations for CBS (default=1000)
--p_seg | p-value for CBS segmentation (default=0.05)
--p_val | p-value for segment validation (default=0.01)
--quantile | Quantile for segment merging (default=0)
--min_ploidy | Minimum ploidy for fitting (default=1.7)
--max_ploidy | Maximum ploidy for fitting (default=3.7)
--min_cellularity | Minimum cellularity/purity for fitting (default=0)
--max_cellularity | Maximum cellularity/purity for fitting (default=1)
--distance_function | Distance function for fitting (RMSD or MAD, default=RMSD)
--threads | Number of threads to use (default=max available)
--mapq | Minimum MAPQ for reads (default=5)
--goi | Path to gene list of interest (BED format)
--size_select | Enable read size selection
--min_read_size | Minimum read sizes for selection
--max_read_size | Maximum read sizes for selection
--no_plot_points | Number of points to plot (default=7500)

## Output Files

**COPYBARA** generates several output files in the specified output directory:

- `{sample}_raw_read_counts.tsv` - Raw read counts per bin
- `{sample}_read_counts_pon_log2r_segmented.tsv` - Segmented and normalised (self, PoN or matched normal) log2 ratios
- `{sample}_fitted_purity_ploidy.tsv` - Fitted purity (tumour fraction) and ploidy values
- `{sample}_segmented_absolute_copy_number.tsv` - Final absolute copy number segments
- `{sample}_copy_number_plot.pdf` - Genome-wide copy number plot

For **COPYBARA focal** analysis the following output files are generated for each gene or region of interest:

- `{sample}_focal_analysis_stats_{roi}.tsv` - Detailed statistics per ROI
- `{sample}_focal_analysis_stats_{roi}.pdf` - COPYBARA focal plot 
In addition, for each sample a summary table (`{sample}_focal_analysis_summary_stats.tsv`) across all ROIs is outputted in the specified output directory.

## Troubleshooting

Please raise a GitHub issue if you encounter issues installing or using COPYBARA.

## License

Apache 2.0 License

Copyright (c) 2024 - All rights reserved.

## Contacts

Carolin Sauer: csauer@ebi.ac.uk

Isidro Cortés Ciriano: icortes@ebi.ac.uk