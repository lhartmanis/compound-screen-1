# Variance adjusted t-test for analyzing differential gene expression
Script used to perform differential expression on data mini-bulk expression data.

The test performed is a variance adjusted t-test between compound treated conditions and negative control conditions. Control samples can be randomized to reduce systemic biases associated with using the same control group for each comparison.

The variance adjustment controls for low replicate numbers in the treatment conditions by first fitting a curve between the mean expression and coefficient of variation for each treatment. Condtions with lower variance than predicted by the regression will be adjusted to the predicted value.

## Dependencies
```
pandas
numpy
scipy
matplotlib
seaborn
sklearn
statsmodels
```

## Usage

`variance_adjusted_ttest.py [-h] -i INDATA_PATH -m META_PATH -o OUTFOLDER -r RNA_TYPE
[-e EQUAL_VARIANCE] 
[-l LABEL_TIME] 
[-c COMPOUND_TIME] [-p PSEUDOCOUNT] [-v VAR_CORRECTION_TYPE] [-t TEST_DICT_PATH]
[--max_data_points MAX_DATA_POINTS]
[--replicate_column REPLICATE_COLUMN] [--subsample_g1 SUBSAMPLE_G1] [--subsample_g2 SUBSAMPLE_G2]
[--conds_to_keep CONDS_TO_KEEP]`


**arguments:**
```
  -h, --help            show this help message and exit
  -i INDATA_PATH, --indata_path INDATA_PATH
                        Path to expression data
  -m META_PATH, --meta_path META_PATH
                        Path to meta data file
  -o OUTFOLDER, --outfolder OUTFOLDER
                        Path to outfolder
  -r RNA_TYPE, --rna_type RNA_TYPE
                        [New/ Total/ Old]
  -e EQUAL_VARIANCE, --equal_variance EQUAL_VARIANCE
                        Assume equal variances in t-test? Default = infer
  -l LABEL_TIME, --label_time LABEL_TIME
                        Label time to analyze. Default = 1hr
  -c COMPOUND_TIME, --compound_time COMPOUND_TIME
                        Compound treatment time to analyze. Default = 1hr
  -p PSEUDOCOUNT, --pseudocount PSEUDOCOUNT
                        What pseudocount to use. Default = 0.001
  -v VAR_CORRECTION_TYPE, --var_correction_type VAR_CORRECTION_TYPE
                        How to correct variance? [both/ linear/ lowess]. Lowess correction is slow with many data points owing to the need to locally fit the regression function. Default = both
  -t TEST_DICT_PATH, --test_dict_path TEST_DICT_PATH
                        Path to dictionary specifying the tests to run.
  --max_data_points MAX_DATA_POINTS
                        How many datapoints to use when training the regressors? More points will make the lowess regression slow. Default = 50000
  --replicate_column REPLICATE_COLUMN
                        What column in the meta file will be used to specify replicates for regression training and DE testing?
  --subsample_g1 SUBSAMPLE_G1
                        Subsample this many conditions from group 1
  --subsample_g2 SUBSAMPLE_G2
                        Subsample this many conditions from group 2
  --conds_to_keep CONDS_TO_KEEP
                        Path to file containing barcodes that pass initial QC

```
## Input

As input, the script uses an expression table, along with a meta data file and specifications of what columns in the meta file to use. An example meta file is provided at `compound-screen/differential_expression/meta.txt`. The script also needs a dictionary specifying the comparisons to be made, an example dictionary in json format can be found at `compound-screen/differential_expression/test_dict.json`.

## Output 

The output folder contains subfolders with differential expression tables for each separate compound used.

## Example

```python3 variance_adjusted_ttest.py -i new_rpkms.txt -o compound_screen/differential_expression/DE_results -r New -v linear -m compound_screen/differential_expression/meta.txt -l 1hr -c 1hr --subsample_g1 10 --test_dict_path compound-screen/differential_expression/test_dict.json --replicate_column compound_name```
