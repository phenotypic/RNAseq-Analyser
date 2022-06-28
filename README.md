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

To install any outstanding libraries, run `R` from the command line, then run `install.packages("package_name")`, replacing `package_name` as necessary.

## Results

After running `individual.R`, you will be presented with the following results:

| Analysis | Description | Figure |
| --- | --- | --- |
| Quality control | Filter cells that have unique feature counts over 2,500 or less than 200, have >5% mitochondrial counts | <img width="580" alt="image" src="https://user-images.githubusercontent.com/33377034/176143563-1065ed69-2acf-407d-b574-3e6c1b43a45e.png"> |
| Identify variable features | Features exhibiting high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others) | <img width="580" alt="image" src="https://user-images.githubusercontent.com/33377034/176143608-c767552a-c1cc-479d-8334-fc136291f754.png">|
| Linear dimension reduction | Perform principal component analysis (PCA) on the scaled data. Useful way to visualise cells and features that define the PCA | <img width="580" alt="image" src="https://user-images.githubusercontent.com/33377034/176143952-8ea82326-b56c-4cb0-84af-82839ca08a1b.png"> |
| Dimensional heat map | Allows for easy exploration of the primary sources of heterogeneity in a dataset | <img width="580" alt="image" src="https://user-images.githubusercontent.com/33377034/176144110-425f32b4-15bd-47ef-9c13-e8185fc995a0.png"> |
| Determine dimensionality | Comparing the distribution of p-values for each PC with a uniform distribution | <img width="580" alt="image" src="https://user-images.githubusercontent.com/33377034/176144299-7e867153-aa77-44f1-b73c-f53125fa3619.png"> |
| Non-linear dimensional reduction | Visualise and explore dataset by grouping similar cells together in low-dimensional space | <img width="580" alt="image" src="https://user-images.githubusercontent.com/33377034/176144340-ec5dd981-a524-4ca6-b028-86625340498b.png">

`compare.R` will output lists of the most highly differentially expressed genes and generate some figures:

| Title | Figure |
| --- | --- |
| Integrated analysis of all cells: clustering | <img width="1695" alt="image" src="https://user-images.githubusercontent.com/33377034/176149122-f6579533-9013-4cc5-ba1f-9e7d8fadf1aa.png"> |
| Visualising spatial information for the most differentially expressed genes | <img width="1305" alt="image" src="https://user-images.githubusercontent.com/33377034/176149310-bb7583ec-ded9-4c35-b367-db3c2674d91c.png"> |
