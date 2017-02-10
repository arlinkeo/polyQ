# Selecting co-expressed genes by threshold in regions associated with HD (Coppen et al. 2016)

setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression/HD_masks_Coppen2016")
library(WGCNA)
library("RDAVIDWebService")
options(stringsAsFactors = FALSE)

#Prepare data and functions
load("../resources/polyQ.RData")
load("meanCor_HDnetworkBD.RData")
probeInfo <- read.csv("../ABA_human_processed/probe_info_2014-11-11.csv")
entrezId2Name <- function (x) { row <- which(probeInfo$entrez_id == x); probeInfo[row, 4]} #Input is single element
make.italic <- function(x) {as.expression(lapply(x, function(x) bquote(italic(.(x)))))}

# Select genes by co-expr. threshold 0.5
selection <- lapply(pQEntrezIDs, function(x){c(x, names(which(meanCor[x,] > 0.50)))})
names(selection) <- pQEntrezIDs
save(selection, file = "../resources/genesets_threshold050_HDregion.RData")

# Functional enrichment of gene sets
david<-DAVIDWebService$new(email="D.L.Keo@tudelft.nl", 
                           url="https://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL", "KEGG_PATHWAY", "BIOCARTA", 
                                 "BIOGRID_INTERACTION", "DIP", "INTACT", "MINT"))
bg_list <- probeInfo$entrez_id # Background list from AHBA probe info
bg <- addList(david, bg_list, idType = "ENTREZ_GENE_ID", listName = "AHBA background", listType = "Background")
bg
setTimeOut(david, 150000)
t <- 0.05 # EASE p-value threshold
info <- lapply(selection, function(geneSet){
  if (length(geneSet) > 1){
    pqname <- entrezId2Name(geneSet[1])
    listname <- paste("HDregions_", pqname, sep = "")
    print(paste("List name:", listname))
    result <- addList(david, geneSet, idType = "ENTREZ_GENE_ID", listName = listname, listType = "Gene")
    print(result)
    setCurrentBackgroundPosition(david, 1)
    fname <- paste("goterms050_", "HDregions_", pqname, ".txt", sep = "")
    getFunctionalAnnotationChartFile(david, fname, threshold=t, count=2L)
    result
  }
})
david
info <- unlist(info, recursive = FALSE)
names(info) <- sapply(names(info), function(s) {arr <- unlist(strsplit(s, "\\.")); paste(arr[1], entrezId2Name(arr[2]), sep = "_")})
info
save(info, file = "../resources/davidgo_geneset_info050_HDregion.RData")