# Load packages
library(affy)
library(hgu133plus2.db)
library(sva)
library(AnnotationDbi)
library(GSEABase)
library(dplyr)
library(tidyverse)
library(qusage)

#Get working directory
setwd("/projectnb/bf528/users/vangogh2022/project_1/biologist/")

#Read differential expression file
sample_data <- read.csv("/projectnb/bf528/users/vangogh2022/project_1/biologist/expression_data.csv", row.names = 1, header = TRUE)

#Order by descending t-value
desc_t <- sample_data %>% arrange(desc(t_score))

#Map probeset IDs to gene symbol with AnnotationDbi
sample_keys <- AnnotationDbi::select(hgu133plus2.db, keys = (row.names(desc_t)), columns = ("SYMBOL"))

#Remove duplicate probeset IDs
remdup <- sample_keys[!duplicated(sample_keys[1]),]

#cbind joins columns in differential expression table and table with mapped probesetID:Symbol
adata <- cbind(remdup, desc_t)

#filters for most significant/smallest p values
newsam <- adata %>%
  group_by(SYMBOL) %>%
  filter(p_adj == min(p_adj)) %>%
  ungroup(SYMBOL)

#Removes null values in Symbol from table
newsam <- newsam[!is.na(newsam$SYMBOL),]

#### DELIVERABLE 1 ####
#take top 1000 and then top 10 up and down regulated genes
top1000_up <- head(newsam, 1000)
top10_up <- head(newsam, 10)
top1000_down <- head(newsam %>% arrange(t_score), 1000)
top10_down <- head(newsam %>% arrange(t_score), 10)
write.csv(top10_up, file = "top10_upregulated_genes.csv", row.names = FALSE)
write.csv(top10_down, file = "top10_downregulated_genes.csv", row.names = FALSE)

#Load in Gene Sets
hallmarks <- getGmt("/usr4/bf528/hmqure/project1_biologist/h.all.v7.5.1.symbols.gmt")
kegg <- getGmt("/usr4/bf528/hmqure/project1_biologist/c2.cp.kegg.v7.5.1.symbols.gmt")
go <- getGmt("/usr4/bf528/hmqure/project1_biologist/c5.go.v7.5.1.symbols.gmt")

#### DELIVERABLE 2: number of genesets in geneset database ####
hallmark_len <- length(names(hallmarks)) #50
go_len <- length(names(go)) #10402
kegg_len <- length(names(kegg)) #186


# Obtain genes which were not expressed
notexp_up <- subset(newsam, !newsam$SYMBOL %in% top1000_up$SYMBOL)
notexp_down <- subset(newsam, !newsam$SYMBOL %in% top1000_down$SYMBOL)

# Defining fisher test function with genelist, geneset, nde= not differentially expressed
fishertest <- function(genelist, geneset, nde)           
{ de_in_gs <- length(intersect(genelist,geneset))    #de_in_gs: differentially expressed genes that are in geneset
de_notin_gs <- length(genelist) - de_in_gs           #de_notin_gs: differentially expressed genes that are not in geneset 
notde_in_gs <- length(intersect(nde,geneset))        #notde_in_gs: not expressed but in geneset
notde_notin_gs <- length(nde) - notde_in_gs          #notde_notin_gs: not differentially expressed and not in geneset
return(c(de_in_gs,de_notin_gs,notde_in_gs,notde_notin_gs))}   #return fishertest values
values

# Initialize data frame that stores results from each gene set after fisher test
kegg_fisher <- data.frame(setname = character(), pval = numeric(), est = numeric(), exp = character(), stringsAsFactors = FALSE)
go_fisher <- data.frame(setname = character(), pval = numeric(), est = numeric(), exp = character(), stringsAsFactors = FALSE)
hallmark_fisher <- data.frame(setname = character(), pval = numeric(), est = numeric(), exp = character(), stringsAsFactors = FALSE)


# Get fisher test values by using for loop, insert into appropriate dataframe for each gmt file

# kegg
for (i in 1:length(kegg))
{
  geneid <- geneIds(kegg[i])
  fisher_up <- fishertest(top1000_up$SYMBOL, geneid[[names(geneid)]], notexp_up$SYMBOL)
  fisher_down <- fishertest(top1000_down$SYMBOL, geneid[[names(geneid)]], notexp_down$SYMBOL)
  upregulated <- fisher.test(matrix(fisher_up,nrow=2))
  downregulated <- fisher.test(matrix(fisher_down, nrow=2))
  kegg_fisher[nrow(kegg_fisher) +1, ] <- c(names(geneid), upregulated$p.value, upregulated$est, 'Up')
  kegg_fisher[nrow(kegg_fisher) +1, ] <- c(names(geneid), downregulated$p.value, downregulated$est, 'Down')}

kegg_fisher <- kegg_fisher %>% mutate(pval = as.numeric(pval), est = as.numeric(est))

# go
for (i in 1:length(go))
{
  geneid <- geneIds(go[i])
  fisher_up <- fishertest(top1000_up$SYMBOL, geneid[[names(geneid)]], notexp_up$SYMBOL)
  fisher_down <- fishertest(top1000_down$SYMBOL, geneid[[names(geneid)]], notexp_down$SYMBOL)
  upregulated <- fisher.test(matrix(fisher_up,nrow=2))
  downregulated <- fisher.test(matrix(fisher_down, nrow=2))
  go_fisher[nrow(go_fisher) +1, ] <- c(names(geneid), upregulated$p.value, upregulated$est, 'Up')
  go_fisher[nrow(go_fisher) +1, ] <- c(names(geneid), downregulated$p.value, downregulated$est, 'Down')}

go_fisher <- go_fisher %>% mutate(pval = as.numeric(pval), est = as.numeric(est))

# hallmarks
for (i in 1:length(hallmarks))
{
  geneid <- geneIds(hallmarks[i])
  fisher_up <- fishertest(top1000_up$SYMBOL, geneid[[names(geneid)]], notexp_up$SYMBOL)
  fisher_down <- fishertest(top1000_down$SYMBOL, geneid[[names(geneid)]], notexp_down$SYMBOL)
  upregulated <- fisher.test(matrix(fisher_up,nrow=2))
  downregulated <- fisher.test(matrix(fisher_down, nrow=2))
  hallmark_fisher[nrow(hallmark_fisher) +1, ] <- c(names(geneid), upregulated$p.value, upregulated$est, 'Up')
  hallmark_fisher[nrow(hallmark_fisher) +1, ] <- c(names(geneid), downregulated$p.value, downregulated$est, 'Down')}

hallmark_fisher <- hallmark_fisher %>% mutate(pval = as.numeric(pval), est = as.numeric(est))

# Using the FDR method to adjust the pvalue in each fisher result
kegg_fisher$FDR <- p.adjust(kegg_fisher$pval, method = "BH", n = length(kegg_fisher$pval))
write.csv(kegg_fisher, "kegg_FDR.csv")

go_fisher$FDR <- p.adjust(go_fisher$pval, method = "BH", n = length(go_fisher$pval))
write.csv(go_fisher, "go_FDR.csv")   

hallmark_fisher$FDR <- p.adjust(hallmark_fisher$pval, method = "BH", n = length(hallmark_fisher$pval))
write.csv(hallmark_fisher, "hallmark_FDR.csv")


#### DELIVERABLE 3 and 4 ####
#Finding significantly enriched genesets 

# kegg

enriched_kegg <- kegg_fisher[kegg_fisher$pval<0.05,]
enriched_kegg_len <- length(enriched_kegg$setname) #21 #counts number of sig enriched gs
enriched_kegg_desc <- enriched_kegg %>% arrange(pval)
enriched_kegg3 <- head(enriched_kegg_desc, 3) #top 3 sig enriched gs
view(enriched_kegg_len)

# go

enriched_go <- go_fisher[go_fisher$pval<0.05,]
enriched_go_len <- length(enriched_go$setname) #888
enriched_go_desc <- enriched_go %>% arrange(pval)
enriched_go3 <- head(enriched_go_desc, 3)
view(enriched_go_len)

# hallmarks

enriched_hallmarks <- hallmark_fisher[hallmark_fisher$pval<0.05,]
enriched_hallmarks_len <- length(enriched_hallmarks$setname) #12 
enriched_hallmarks_desc <- enriched_hallmarks %>% arrange(pval)
enriched_hallmarks3 <- head(enriched_hallmarks_desc, 3) 
view(enriched_hallmarks_len)


write.csv(enriched_kegg3, file="top3_enriched_kegg.csv")
write.csv(enriched_go3, file="top3_enriched_go.csv")
write.csv(enriched_hallmarks3, file="top3_enriched_hallmarks.csv")




