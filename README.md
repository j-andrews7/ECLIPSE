<p align="center" width="100%">
    <img src="inst/logo/ECLIPSE_Hex.png" alt="ECLIPSE" height="330">
</p>

---

ECLIPSE (**E**nhancer **C**alling and **L**inking with **I**ntegrated **P**rofiling and **S**tructure **E**valuation) provides a performant 
implementation of the [rank ordering of super enhancers (ROSE)](http://younglab.wi.mit.edu/super_enhancer_code.html) method for identifying super enhancers.
It provides options to increase the robustness of ROSE via signal transformations prior to thresholding and additional thresholding approaches.
It also increases flexibility by exposing parameters hidden in the original implementation.
ECLIPSE additionally contains novel functionality to identify sub-structural changes within enhancer regions between groups via sliding window and binning approaches.
It also contains visualization functions to generate high-quality plots of specific loci alongside arbitrary user-provided data.

This project was conceptualized for and initially developed during the [St. Jude Children's Research Hospital KIDS24 BioHackathon](https://github.com/stjude-biohackathon) by:
- Jared Andrews (team lead)
- Alyssa Obermayer
- Nicolas Peterson
- Kelly McCastlain
- Jacqueline Myers
- Avery Bradley

It even snagged a prize for "Most Technically Impressive" project.

**This package is under active development and may break at any time.**

## Installation

Currently, the package can be installed from GitHub:

```
library(devtools)
devtools::install_github("stjude-biohackathon/ECLIPSE")
```

## Usage

Given paths to a BAM file and a BED file of peaks, ROSE can be run with the `run_rose` function.
Optionally, a control BAM file for input or IgG from the same sample can be provided.

Alternatively, `bamFile` objects for the treatment and control signal and a `GRanges` object for the peaks can be provided.

The output is a `GRanges` object containing all putative enhancers with their super enhancer designation in the `super` metadata column.

Below is an example of running ROSE on a BAM file of H3K27ac MINT ChIP-seq, an [input control](https://www.encodeproject.org/experiments/ENCSR056PPJ/) BAM file, and a BED file of peaks from [this ENCODE experiment of human naive B cells](https://www.encodeproject.org/experiments/ENCSR660EVU/).

```r
# We'll use the BiocFileCache package to download and cache the files, which will take a few minutes the first time they're used.
library(BiocFileCache)
bcf <- BiocFileCache(ask = FALSE)
treat_url <- "https://www.encodeproject.org/files/ENCFF993DJI/@@download/ENCFF993DJI.bam"
treat_path <- bfcrpath(bfc, treat_url)

control_url <- "https://www.encodeproject.org/files/ENCFF821MAI/@@download/ENCFF821MAI.bam"
control_path <- bfcrpath(bfc, control_url)

peaks_url <- "https://www.encodeproject.org/files/ENCFF590DFY/@@download/ENCFF590DFY.bed.gz"
peaks_path <- bfcrpath(bfc, peaks_url)

naiveB1_enhancers <- run_rose(treatment = treat_path, control = control_path, peaks = peaks_path)

naiveB1_enhancers
```

## Development Roadmap

- ~~Add missing ROSE functionality.~~
  - ~~TSS exclusion from stitching process.~~
  - ~~Overlap of TSS of 3 or more unique genes canceling stiching for a putative enhancer. ~~
    - ~~And ability to disable this process.~~
  - ~~Enhancer-gene annotations (within 50 kb by default for ROSE).~~
    - ~~Ability to limit to expressed genes.~~
- Allow BigWig or GRanges signal as input.

## References

If you use this package in published research, please cite the original papers utilizing and describing ROSE:

[Master Transcription Factors and Mediator Establish Super-Enhancers at Key Cell Identity Genes](http://www.cell.com/abstract/S0092-8674(13)00392-9)
Warren A. Whyte, David A. Orlando, Denes Hnisz, Brian J. Abraham, Charles Y. Lin, Michael H. Kagey, Peter B. Rahl, Tong Ihn Lee and Richard A. Young
Cell 153, 307-319, April 11, 2013

[Selective Inhibition of Tumor Oncogenes by Disruption of Super-enhancers](http://www.cell.com/abstract/S0092-8674(13)00393-0)
Jakob Lovén, Heather A. Hoke, Charles Y. Lin, Ashley Lau, David A. Orlando, Christopher R. Vakoc, James E. Bradner, Tong Ihn Lee, and Richard A. Young
Cell 153, 320-334, April 11, 2013
