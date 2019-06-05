# Script for downloading data used in 
# Blomstedt et al.: "Modelling-based experiment retrieval: 
# A case study with gene expression clustering".
# Author(s): Sohan Seth, Paul Blomstedt 


######################################################
#### Create folder structure #########################
currentDir <- system2("pwd", stdout=TRUE, stderr=TRUE)
setwd(currentDir)
dir.create(path = 'Data')
dir.create(path = 'Data/TAB')
dir.create(path = 'Data/TAB/Raw')
dir.create(path = 'Data/TAB/Analysed')
dir.create(path = 'Data/EFO')
dir.create(path = 'Data/MAT')


###############################################################################
#### Get list of accession numbers for AFFY-44 experiments ####################

library("SPARQL")

### Get accessions
# Extract experiments, description and respective platforms (AffyMetrix)
endpoint <- "http://www.ebi.ac.uk/rdf/services/atlas/sparql"
query <- "PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX owl: <http://www.w3.org/2002/07/owl#>
PREFIX dcterms: <http://purl.org/dc/terms/>
PREFIX obo: <http://purl.obolibrary.org/obo/>
PREFIX sio: <http://semanticscience.org/resource/>
PREFIX efo: <http://www.ebi.ac.uk/efo/>
PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/>
PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/>
SELECT DISTINCT ?experiment ?accession ?description ?platform WHERE
{?experiment
a atlasterms:Experiment ;
dcterms:identifier ?accession ;
dcterms:description ?description ;
atlasterms:hasAnalysis
[atlasterms:hasExpressionValue 
[atlasterms:isMeasurementOf 
[atlasterms:partOfPlatform ?platform
]
]
]
filter regex (?platform, \"A-AFFY-44\", \"i\")
}"
data <- SPARQL(endpoint,query)


############################################################################
################# Download gene expression data ############################
setwd('Data/TAB/Raw')
for (accession in data$results$accession ){
  system(paste('curl -o ',accession,'.zip \"http://www-test.ebi.ac.uk/gxa/directDownload/experimentExpressions?eacc=',accession,'\" && unzip ',accession,'.zip',sep=''))
}
setwd(currentDir)

############################################################################
################# Download analysed data ###################################
setwd('Data/TAB/Analysed')
for (accession in data$results$accession ){
  system(paste('curl -o ',accession,'.zip \"http://www-test.ebi.ac.uk/gxa/directDownload/experimentAnalytics?eacc=',accession,'\" && unzip ',accession,'.zip',sep=''))
}
setwd(currentDir)

#####################################################################################
############ Download and extract selected EFOs #####################################
d <- list()
for (accession in data$results$accession)
{
  query <- paste("PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
                 PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
                 PREFIX owl: <http://www.w3.org/2002/07/owl#>
                 PREFIX dcterms: <http://purl.org/dc/terms/>
                 PREFIX obo: <http://purl.obolibrary.org/obo/>
                 PREFIX sio: <http://semanticscience.org/resource/>
                 PREFIX efo: <http://www.ebi.ac.uk/efo/>
                 PREFIX atlas: <http://rdf.ebi.ac.uk/resource/atlas/>
                 PREFIX atlasterms: <http://rdf.ebi.ac.uk/terms/atlas/>
                 SELECT DISTINCT ?experiment ?description ?EFOType ?EFOValue WHERE
{?experiment
                 a atlasterms:Experiment ;
                 dcterms:identifier ?accession ;
                 dcterms:description ?description ;
                 atlasterms:hasAssay
                 [atlasterms:hasSample
                 [atlasterms:hasSampleCharacteristic
                 [atlasterms:propertyType ?EFOType ;
                 atlasterms:propertyValue ?EFOValue
                 ]
                 ]
                 ]
                 filter regex (?accession, \"",accession,"\", \"i\")
}", sep = "")
  tmp <- SPARQL(endpoint,query)
  for (cat in names(tmp$result)){
    d[[cat]] <- c(d[[cat]], tmp$result[[cat]])
  }
}

#### transform into data frame and save as a csv file
EFOdata <- data.frame(matrix(unlist(d), nrow=length(d[[1]]), byrow=FALSE))
colnames(EFOdata) <- attributes(d)$names
write.csv(EFOdata, file = 'Data/EFO/EFOdata.csv',row.names=FALSE)

#### extract information on ORGANISM PART and save as csv file
organism_part <- subset(EFOdata,EFOType == "organism part", select = c(experiment,EFOValue))
organism_part <- droplevels(organism_part)
organism_part$experiment <- sub(".*?atlas/(.*?)>.*", "\\1", organism_part$experiment) # extract experiment identifier
write.csv(organism_part, file = 'Data/EFO/organism_part.csv',row.names=FALSE)

#### extract information on CELL TYPE and save as csv file
cell_type <- subset(EFOdata,EFOType == "cell type", select = c(experiment,EFOValue))
cell_type <- droplevels(cell_type)
cell_type$experiment <- sub(".*?atlas/(.*?)>.*", "\\1", cell_type$experiment) # extract experiment identifier
write.csv(cell_type, file = 'Data/EFO/cell_type.csv',row.names=FALSE)

#### extract information on DISEASE and save as csv file
disease <- subset(EFOdata,EFOType == "disease", select = c(experiment,EFOValue))
disease <- droplevels(disease)
disease$experiment <- sub(".*?atlas/(.*?)>.*", "\\1", disease$experiment) # extract experiment identifier
write.csv(disease, file = 'Data/EFO/disease.csv',row.names=FALSE)



