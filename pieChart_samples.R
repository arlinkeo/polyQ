# pie chart of number of samples
##### Obtain binary vectors indicating presence of sample across regions and donors #####
setwd("C:/Users/dkeo/surfdrive/polyQ_coexpression")
options(stringsAsFactors = FALSE)

load("resources/sampleIDs.RData")
sampleIDs["cerebellum"] <- NULL
ontology <- read.csv("../ABA_human_processed/Ontology_edited.csv")

nSamples <- sapply(sampleIDs, function(r){
  sapply(r, sum)
})

structures <- names(sampleIDs)[-c(which(names(sampleIDs) == "HD_region"))]
structures <- sapply(structures, function(s)gsub("_", " ", s))
colorOther <- "grey"
colorsStruct <- sapply(ontology$color_hex_triplet[match(structures, ontology$name)], function(c){paste0("#", c)})
colorsStruct <- c(colorsStruct[-1], colorOther, colorsStruct[1])
colorsHD <- c("#956b53", colorOther)

sumDonors <- apply(nSamples, 2, sum)
nonHD <- sumDonors["brain"]-sumDonors["HD_region"]
vecHD <- c(sumDonors["HD_region"], nonHD)
names(vecHD) <- c("HD-associated region", "rest")
names(vecHD) <- paste0(names(vecHD), " (", vecHD, ")")
sum(vecHD)

sumStruct <- sumDonors[-which(names(sumDonors) %in% c("brain", "HD_region"))]
nonStruct <- sumDonors["brain"] - sum(sumStruct)
vecStruct <- c(sumStruct, nonStruct)
names(vecStruct) <- c(structures[-1], "rest")
names(vecStruct) <- paste0(names(vecStruct), " (", vecStruct, ")")

pdf(file = "Images/pieChart_samples.pdf", 12, 5)
layout(matrix(c(1:2), 1, 2))
par(mai = c(0,1,0,1), oma = c(0,2,0,2))
pie(vecHD, clockwise = TRUE, col = colorsHD, radius = 0.5, main = "Total number of samples")
pie(vecStruct, clockwise = TRUE, col = colorsStruct, radius = 0.5)
dev.off()