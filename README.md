# GSR: a manually curated and optimised taxonomical database for 16S rRNA amplicon analysis

<img src="GSR_logo.png">

## What is the GSR database?

GSR database (Greengenes, SILVA, and RDP database) is an integrated and manually curated database for bacterial and archaeal 16S amplicon taxonomy analysis. Unlike previous integration approaches, this database creation pipeline includes a taxonomy unification step to ensure consistency in taxonomical annotations. The database was validated with three mock communities and two real datasets and compared with existing 16S databases such as Greengenes, GTDB, ITGDB, SILVA, RDP, and MetaSquare. Results showed that the GSR database enhances taxonomical annotations of 16S sequences, outperforming current 16S databases at the species level. The GSR database is available for full-length 16S sequences and the most commonly used hypervariable regions: V4, V1-V3, V3-V4, and V3-V5. 

## Usage guidelines

### Recommended usage
1. Download the GSR database of the region of interest from [here](https://manichanh.vhir.org/gsrdb/).
2. Train a QIIME2 classifier with the database and use it to get the taxonomy profile (see [QIIME2 tutorial](https://docs.qiime2.org/2023.2/tutorials/moving-pictures-usage/) for more information). 

<!-- we recommend following the QIIME2 pipeline to preprocess the reads??? -->

### Example

In this example we will obtain the taxonomy profile of our 16S V4 amplicon reads. Example files: `blablabla.fastq etc` 

1. Install QIIME2 <!-- for command line or for python????-->
2. Read preprocessing 
	1. Demultiplex
	```

	```
	2. Denoise
	```

	```

3. Download GSR V4 database
```
wget database download link
```

4. Train naive-bayes classifier with default parameters
```

```
5. Get taxonomy profile, setting confidence level to 'disable'. 
```

```
<!-- explain 
	1. how to download  
	2. how to use it with qiime2 (explain also qiime2 preprocessing of raw reads to get repseqs) (OR refer to qiime2 tutorial) (OR both!)
	IMPORTANT: say the recommended parameters for qiime classifier-->
<!-- upload also the trained classifiers???? or let people train them ???? 
	trained classifier for full16S database is too big for upload :(( 
	maybe we can keep the database downloading in the webpage!! -->

## Citing GSR database 

If you use GSR database, please cite:

*Leidy Alejandra G. Molano, Sara Vega-Abellaneda, Chaysavanh Manichanh. GSR: a manually curated and optimised taxonomical database for 16S rRNA amplicon analysis. UNPUBLISHED YET*