# RNAseq-Analyser

The scripts in this repository are used to analyse RNA-seq data using [Seurat](https://github.com/satijalab/seurat). `individual.R` allows you to analyse a single dataset, while `compare.R` lets you compare two datasets to uncover differential gene expression.

Both scripts expect dataset folders to include the following correctly formatted files: `barcodes.tsv`, `genes.tsv`, `matrix.mtx`. Take a look at the [RNAseq-Parser](https://github.com/phenotypic/RNAseq-Parser) repository if you need to format your input files.

The sample data referred to in this repository can be obtained from [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192498). It must be processed using the RNAseq-Parser repository before use.

## Usage

Download with:
```
git clone https://github.com/phenotypic/RNAseq-Analyser.git
```

To run the scripts you must have [RStudio](https://www.rstudio.com/) installed. Download through [brew](https://brew.sh/) by running: `brew install --cask rstudio`

Here are the plots you can expect from `individual.R`:


Here are the plots you can expect from `compare.R`:
