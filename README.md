# RExGC
Retrieval of experiments using gene clustering (RExGC)

This file describes the workflow for downloading, processing and analysing the data in:<br/> 
<p>
    Blomstedt, P., Dutta, R., Seth, S., Brazma, A. and Kaski, S. (2016).
    <b>Modelling-based experiment retrieval: A case study with gene expression clustering.</b>
    <i>Bioinformatics</i>, <b>32</b>(9), 1388&ndash;1394.   
    <a href="http://bioinformatics.oxfordjournals.org/content/32/9/1388">[Link]</a>
</p>

The files in this repository are also available as a [zip file](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/32/9/10.1093_bioinformatics_btv762/4/btv762_Supplementary_Data.zip?Expires=1560340402&Signature=TG~QUE7J8hVpkbX5wcWqvwKcXNYlWg9yKLoznl5NPjKMP4XWDKIHIyJFRMzjEspZ93mBaNw8oz5UfQajjyDkG-CoYPExGafJylpbPG9vlgg8YTGF-D8l6XqSLLZj7U~tbigZitus-3JPUvGMqkRwvAg4QVn1t-tKdpqx67k-ljtJqejMikKumA13OFmJUoNP7oFNHEAL13n33g00vzvPz8kPJon3d0sJ-hc9orqp64k9pI9QJ~AfHekv772xu07V8uqmWttj1Zc46w7Men7dAvsSQMPsQFx2Cyqu3AWrz2Zqmb9ohZy0tpU0K98IL74t1t4zNukU6cQvWfj349BFkQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA) from the publisher's webpage. 
The code runs in MATLAB (version R2012b or later) and has been designed and tested in linux.
For a demo, run the `demo.m` script. 
Datasets are stored in the Data folder with the following subfolders 
* EFO: csv files of EFO data 
* MAT: mat files generated as output from the MATLAB functions below
* TAB: tab files of gene expression data downloaded from Expression Atlas. The tab files (~20GB) are not prestored 
	due to space considerations but can be downloaded using the R script provided in step 1. 
 
NB! Step 1 is a preliminary step for downloading and preprocessing the data and may take a considerable amount of 
time to run. Steps 2 and 4b require the availability of the data downloaded in step 1. To use preprocessed datasets
based on downloads from 4th Jun 2014, steps 1, 2 and 4b may be omitted.


## 1. Create folder structure, download and preprocess data 

**Function**:<br/>
`createData;`

**Output**:<br/>
- EFO data: `EFOdata.csv`, `cell_type.csv`, `disease.csv`, `organism_part.csv`
- Ground truth matrices: `G44_EFO_combined.mat`, `G44_cell_type.mat`, `G44_disease.mat`, `G44_organism_part.mat`; 
	these files also contain lists of accession numbers for experiments used in the analyses as well as 
the annotations used to construct the gold standard
- tab-files of measurement and analytics data for gene expression experiments indentified by accession number

**Requires**:<br/>
Downloading and initial preprocessing of the data is done in R and requires the
SPARQL package to be installed.


## 2. Selection of genes

**Function**:<br/>
`selectGenes(nTop);`

**Input**:<br/>
`nTop` - number of genes to include for each experiment from the top of the default ordering given in the tab files 
of analytics data 

**Output**:<br/>
`['44_top_', num2str(nTop) '.mat']` - processed datasets for analysis, also contains a list of selected genes and 
the full set of all genes

**Requires**:<br/>
Data in tab files and `G44_EFO_combined.mat`

**Example**:<br/>
`selectGenes(5);`


## 3. Clustering the data

**Function**:<br/>
`rexClustering(algorithm_type,nTop);`

**Input**:<br/>
- `algorithm_type` - type of clustering algorithm used: `'PPM'` (full model space), `'PPMkm'` (reduced model space) or 
	`'k-means'` (fixed no. of clusters)
- `nTop` - identifier for the dataset corresponding to `nTop` (as explained above)

**Output**:<br/>
`[algorithm_type '_clust_top_' num2str(nTop) '.mat']` - vector of cluster labels

**Requires**:<br/>
`['44_top_', num2str(nTop) '.mat']` 

**Example**:<br/>
`rexClustering('PPM',5);`


## 4a. Evaluate similarity between obtained clusterings using Normalized Information Distance (NID)

**Function**:<br/>
`computeNid(algorithm_type,nTop);`

**Input**:<br/>
`algorithm_type` and `nTop` specified as above

**Output**:<br/>
`[algorithm_type '_NID_top_' num2str(nTop) '.mat']` - matrix of NIDs between the experiments

**Requires**:<br/>
`[ algorithm_type '_clust_top_' num2str(nTop) '.mat']` 

**Example**:<br/>
`computeNid('PPM',5);`


## 4b. Similarity matrices for baseline approaches (may be omitted if prestored datasets are used)

### Likelihood-based approach

**Function**:<br/>
`evaluateLogML;`

**Output**:<br/>
`LogML.mat` - similarity matrix for likelihood-based approach

**Requires**:<br/>
`PPM_clust_top_5.mat`, `44_top_5.mat`


### Differential expression based approach

**Function**:<br/>
`DEpValCorrelation;`

**Output**:<br/>
`Data/MAT/DE.mat` - similarity matrix for likelihood-based approach for the DE approach
(Pearson's correlation)

**Requires**:<br/>
Data in tab files and `G44_EFO_combined.mat`


## 5. Visualize results

**Function**:<br/>
`rexVisualization(EFOtype,algorithm,nTop);`

**Input**:<br/>
- `EFOtype` - EFO type to display results for: `'cell_type'`, `'disease'` or `'organism_part'`
- `algorithm` - cell array of retrieval methods, see example below
- `nTop` - specified as in the steps above

**Output**:<br/>
precision-recall curve
 
**Requires**:<br/>
`['Data/MAT/G44_' EFOtype '.mat']` and `['Data/MAT/' algorithm{i} '_NID_top_' num2str(nTop) '.mat']`

**Example**:<br/>
```
algorithm = {'PPM','LogML','DE'};
rexVisualization('cell_type',algorithm,5);
legStr = [{'Random'},algorithm];
legend(legStr,'Location','best');
```
