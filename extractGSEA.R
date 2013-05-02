
# Extract genes included in pathways

library(GSEABase)
source("/Users/fredcommo/Documents/Stats/Fred_Scripts/GeneRequest.v6.R")


# Extraire les GeneSets
# syntaxe = system.file(“subfolder”, ”file_name.xml”, package=”package_name”)

address <- system.file("datasets", "msigdb_v3.0.xml", package="GSEABase")
BSets<- getBroadSets(address)

p<-length(BSets)
BNames<-data.frame(index=seq(1,p), setname=NA, Cat=NA, SubCat=NA, Nb.Genes=0)
for (i in 1:p){
BNames[i,2]<-BSets[[i]]@setName
BNames[i,3]<-BSets[[i]]@collectionType@category
BNames[i,4]<-BSets[[i]]@collectionType@subCategory
BNames[i,5]<-length(BSets[[i]]@geneIds)
}
# BNames contient les GeneSets, catégorie (voir GSEA) et nombre de gènes
# Afficher BNames[1:10,] pour exemple

table(BNames$Cat)
table(BNames$SubCat)

ExtractGene <-function(GeneSet.list, BSets){
	index<- GeneSet.list$index
	output<-data.frame()
	for (i in 1:length(index)){
		tmp.genes<-	as.vector(BSets[[index[i]]]@geneIds)
		tmp.setname<- as.vector(BSets[[index[i]]]@setName)
		tmp<- cbind.data.frame(setname=tmp.setname, gene= tmp.genes)
		output<-rbind.data.frame(output,tmp)
		}
	output
}

CP <- BNames[which(BNames$SubCat == "CP"), ]
head(CP)

myList <- ExtractGene(CP, BSets)
head(myList)

# Use GeneRequest to get GeneIds
myIds <- GeneRequest.v6(as.character(myList$gene), "gene")
head(myIds)
myGenes <- cbind.data.frame(myIds, Pathway = myList$setname)
head(myGenes)


