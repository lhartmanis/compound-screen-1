# New RNA inference

This folder contains a computational pipeline that infers the fraction of new RNA per gene within a sample.

Describe how to design the yaml and run the pipeline.
Add a note on runtime and memory usage.


The inference script runs in sequential batches to keep memory usage low. The number of genes to be inferred per batch can be set in the Snakemake configuration file. Processing data using smaller batches requires less memory usage and does not significantly alter the runtime.

The pipeline has been benchmarked on a desktop computer with the following specs:
```
Hardware:
CPU: AMD Ryzen Threadripper 3960X, 24 core 48 threads
GPU: NVIDIA GeForce GTX 1080 Ti
Memory: 256 GB 

Drivers:
Driver Version: 470.57.02
CUDA Version: 11.4 
```

Processing data from a large-scale experiment with 10 million data points, running on 30 parallel processes and performing 10.000 MCMC sinulations per batch finished in 96 hours and each batch required 130 GB RAM. 

## Dependencies
```
joblib=v1.0.1
pysam=0.17.0
gtfparse=v1.2.1
pandas=1.1.4 
numpy=1.19.5 
scipy=1.6.3
tensorflow=2.6.0
tensorflow_probability=0.14.0
```

## Usage
Parameters for each step is entered into the `config.yaml` file. An example file can be found at `compound-screen/new_RNA_inference/config.yaml`. Once the yaml file is complete, place it in the desired run directory together with the Snakefile found at `compound-screen/new_RNA_inference/config.yaml`.

The full pipeline can then easily be started through: 

`snakemake --cores 10`

## Pipeline details
The pipeline consists of four separate steps:

### 1) Conversion counting from sequencing data
This step counts the number of observed mismatches to the reference genome observed in each sequencing read and computes probabilities of errors and specific conversions. This is performed by the script `stranded_conversion_rates_and_probs.py`

#### Input
Configuration file in yaml-format.
Conversion counting starts from an aligned BAM-file. An example of how to generate one from fastq files is found in `path_to_file`.

The supplied BAM file needs to contain read tags for genes and condition barcodes. By default, the pipeline looks for the gene tags `GE` and `GI` and barcode tag `BC`.

#### Output
`details.txt` Number of observed mismamtches and observed T:s per condition and gene.

`pc_pe.txt` Computed mismatch probabilities for each condition.

`conversionRates.csv` Conversion rates for each base conversion per condition.

### 2) Combination and indexing of conversion count files
This step combines the output from _step 1_ into the file that will be used as input to step _3_. The file gets indexed with [tabix](http://www.htslib.org/doc/tabix.html) to allow fast data retreival. This step is carried out by the script `prepare_and_index_indata.py`

#### Input
Path to directory where the output from the previous step is stored. By default `./data`

#### Output
```
combined_results.txt.gz
combined_results.txt.gz.tbi
total_comparisons.txt
```

### 3) Inference of $\pi_g$ - the fraction new RNA per gene in a treatment
This is where the rubber hits the road. Starting from the `combined_results` file generated previously, this step infers the fraction of new RNA per gene and condition using an MCMC sampling algorithm. To reduce the amount of RAM needed by the algorithm, it is run in sequential batches processing a set number of genes per iteration before restarting with the next data chunk. The task of handling the automated restart is performed by the script `runner_script.py`, and the inference is performed by the script `pi_g_inference.py`.

This is the slowest part of the pipeline

#### Input
`combined_results.txt.gz`

#### Output
Individual text files are generated for each inference (treatnent-gene pair). The files contain the most likely inferred $\pi_g$ value as well as the 95% credible interval of the inference along with the number of reads used in the inference.

### 4) Collection of outdata and removal of intermediary folders
This step combines the data from each file in _step 3_ into a combined output dataframe

#### Input
Results folders from _step 3_

#### Output
`pi_g_df.txt.gz`

The resulting $\pi_g$ matrix contains the fraction of new RNA inferred per gene in each condition. Multiplying this with an expression matrix generates new RNA expression levels, as shown in this [example notebook](https://github.com/lhartmanis/compound-screen/blob/main/examples/perform_qc_and_separate_expression_matrices.ipynb).

