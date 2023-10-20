# Transcription factor binding enrichment 

Contains scripts and data for analyzing transcription factor binding enrichment in a target gene set compared to background.
Data files containing ChIP-seq data needs to be unzipped before the pipeline can be run.

This script is designed to be used for large-scale compound-screening experiments and requires a determined - although not very complicated - directory structure to work. The top directory, further referred to as `infolder`, contains subdirectories for each treatment. Each subdirectory in turn contains text files for the target and background gene set where one gene name is written per row.

## Usage

`threaded_enrichment_target_vs_bckg.py [-h] [-i INFOLDER] [-p PROCESSES] [--target_name TARGET_NAME] [--background_name BACKGROUND_NAME] [-c CHIP_TYPE] [-s SIG_TEST]`

**arguments**
```
 -h, --help            show this help message and exit
  -i INFOLDER, --infolder INFOLDER
                        Path to a folder containing GSEA data in one or more subfolders. It is a mess and I need to change it soon.
  -p PROCESSES, --processes PROCESSES
                        Number of concurrent processes to utilize
  --target_name TARGET_NAME
                        Name of target gene list. Default is target.txt
  --background_name BACKGROUND_NAME
                        Name of background gene list. Default is background.txt
  -c CHIP_TYPE, --chip_type CHIP_TYPE
                        Decide what ChIP-data to use. Implemented files are: [filtered_peaks]. Default is filtered_peaks
  -s SIG_TEST, --sig_test SIG_TEST
                        Decide which significance test to use. Implemented tests are: [Ksamp_AndersonDarling/ Wilcoxon_ranksums]. Default is Wilcoxon ranksums
```

## Input
Starting from `infolder`, each treatment will have its own subdirectory wherein two gene lists are saved as text files with gene names on each row.

## Output
Output folders contain enrichment of transcription factor binding in each of eight genomic regions for each treatment

```
upstream
upstreamproximal
downstreamproximal
promoter
genebody
polyAsite
proximal
downstreamctrl
```

## Example
`python3 threaded_enrichment_target_vs_bckg.py -i infolder -p 30 --target_name upregulated.txt --background_name downregulated.txt -c filtered_peaks -s Wilcoxon_ranksums`
