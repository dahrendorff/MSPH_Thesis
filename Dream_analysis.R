#1. Standard RNA-seq processing

library('variancePartition')
library('edgeR')
library('BiocParallel')
library('tibble')
library(compositions)
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)


#countmatrix is our featurecounts file named "cleaned_combined_countmatrix_remove_dub.txt"

countMatrix<-read.table("C:\\Users\\jdahr\\Desktop\\dream\\cleaned_combined_countmatrix _remove_dub.csv",header = TRUE,
                        fill = TRUE, check.names=FALSE, sep=",")



#check.names=FALSE gets rid of the X infront of the samplenames. 
#X is added in R when column names are not in the correct form or start with a number or special character.
#make the  entries of column1 as the rownames because the columns can only contain numerical entries.

# convert column to rownames
countMatrix <- column_to_rownames(countMatrix, var="ENSEMBLEID")

head(countMatrix[,1:5])

#the countMatrix is now in the desired format.
#read counts = number of reads that map (i.e align) to each gene
#filter genes by number of counts
#Genes with very low counts across all libraries provide little evidence for differential expression
#These genes should be filtered out prior to further analysis.

isexpr = rowSums(cpm(countMatrix)>0.1) >= 5
View(isexpr)

# Standard usage of limma/voom
#DGEList is an object to hold our read counts. Container for counts themselves sometimes also metadata.

geneExpr = DGEList( countMatrix[isexpr,] )
geneExpr = calcNormFactors( geneExpr )

#Calculate normalization factors is used to scale the raw library sizes
#calcNormFactors doesn't normalize the data, it just calculates normalization factors for use downstream
#you have to use the name of the object that contains multiple lists and$ to specify which list you want to look at.

View(geneExpr$counts)
View(geneExpr$samples)

#metadata option #1
#"metadata" used in the dream vigenette is our"dnhs_fulldata_fordeseq_removed_proportions_with_cellpro_csv" file. It contains the RNAIDs, RespIDs,wave#info
#the cellproportions ect.# cell proportions computationally derived with Cibersortx & added to metadata.

#metadata with actual cellproportions. Cellproportions were  previously merged into the metadatafile.
# 5 fractions sum to 1.0 for each sample means that the fractions really span only 4 dimensions
#(i.e., knowing three determines the fourth) problematic in regression model
metadata<-read.csv("C:\\Users\\jdahr\\Desktop\\dream\\dnhs_dream_metadata_input.csv")
head(metadata)
table(metadata$Sex)
table(metadata$ptsd)
table(metadata$CD4T)

# check data order
table(colnames(countMatrix) == metadata$rnaid)

#Skip this step if you already have merged Metadata file with transformed cell proportions
#go to line #83
#metadata option #2 #Metadata needs to be merged with cell proportions file.
METADATA <-read.csv("C:\\Users\\jdahr\\Desktop\\practice\\dnhs_dream_metadata_input.csv")
df_cellFractions = read.csv("C:\\Users\\jdahr\\Desktop\\practice\\Cell_proportions_NSCLC_PBMC_both_timepoints.csv")
df_cellFractions <- column_to_rownames(df_cellFractions, var="rnaid")

# Transform the 5 cell types into 4 variables using the isometric log ratio "compositions" package ilr function
celFrac_ilr = ilr(df_cellFractions)

colnames(celFrac_ilr) = paste0("cellFrac_ilr_", 1:4)

METADATA1 = merge(METADATA, df_cellFractions, by.x = 'rnaid', by.y='row.names')
METADATA1 = merge(METADATA1, celFrac_ilr, by.x = 'rnaid', by.y='row.names')

write.csv(METADATA1,"C:\\Users\\jdahr\\Desktop\\practice\\METADATA_log_transf_prop.csv")

#figure out if you should delete non-log transf cell proportions from columns.

#shortcut to metadata file
Metadata1 = read.csv("C:\\Users\\jdahr\\Desktop\\dream\\METADATA_log_transf_prop_both.csv")

#you can skip to Dream analysis to Line #124
#2. Limma Analysis
#Limma has a built-in approach for analyzing repeated measures data using duplicateCorrelation(). 
#The model can handle a single random effect, and forces the magnitude of the random effect to be the same across all genes.
# apply duplicateCorrelation is two rounds

design = model.matrix( ~ ptsd, Metadata1)


vobj_tmp = voom( geneExpr, design, plot=TRUE)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=Metadata1$resp)




# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
# Otherwise, use the results from the first voom run

vobj = voom( geneExpr, design, plot=FALSE, block=Metadata1$resp, correlation=dupcor$consensus)
View(vobj_tmp$targets)
View(vobj_tmp$E)

# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 

dupcor <- duplicateCorrelation(vobj, design, block=Metadata1$resp)

# But this step uses only the genome-wide average for the random effect

fitDupCor <- lmFit(vobj, design, block=Metadata1$resp, correlation=dupcor$consensus)

# Fit Empirical Bayes for moderated t-statistics

fitDupCor <- eBayes( fitDupCor )

#3.Dream Analysis
#The dream method replaces two core functions of limma with a linear mixed model.
###Otherwise dream uses the same workflow as limma with topTable(), since any statistical differences are handled behind the scenes.

# Specify parallel processing parameters
# this is used implicitly by dream() to run in parallel

param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)

# The variable to be tested must be a fixed effect

#1.
#first testing cases versus controls with a formula like this
#test if the disease coefficient is different from zero to see if there is a disease effect after removing the variation across waves
#and accounting for the repeated measurement of donors.

form <- ~ ptsd + (1|resp) + wave + Sex + age + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3 + cellFrac_ilr_4

#form without wave
form <- ~ ptsd + (1|resp) + Sex + age + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3 + cellFrac_ilr_4


# estimate weights using linear mixed model of dream

vobjDream = voomWithDreamWeights( geneExpr, form, Metadata1 )

# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test

fitmm = dream( vobjDream, form, Metadata1 )


#fitmm = dream( vobjDream, form, Metadata1, ddf="Kenward-Roger" )

# Examine design matrix

head(fitmm$design, 3)

# Get results of hypothesis test on coefficients of interest

topTable( fitmm, coef='ptsd', number=3 )





#safe fitmm data
save(fitmm, file = "C:\\Users\\jdahr\\Desktop\\dream\\fitmm.Rdata")


topTable( fitmm, number=5 )

#What genes are most differentially expressed?

top.table <- topTable( fitmm, coef='ptsd', sort.by = "P", number=Inf ) 
top.table <- topTable( fitmm, coef='SexMale', sort.by = "P", number=Inf ) 

#How many DE genes are there?

length(which(top.table$adj.P.Val < 0.05))

#Write top.table to a file

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]
write.table(top.table, file = "C:\\Users\\jdahr\\Desktop\\dream\\TopTableOutput_ptsd_cases_vs_control_nowave.txt", row.names = F, sep = "\t", quote = F)

#2.
#second, test if there is a wave-by-disease interaction using formula
#test the coefficient disease:wave to test of the association between disease
#and gene expression changes between the two waves

form <- ~ ptsd + wave + ptsd:wave + (1|resp) + Sex + age + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3 + cellFrac_ilr_4

# estimate weights using linear mixed model of dream

vobjDream = voomWithDreamWeights( geneExpr, form, Metadata1 )
# Fit the dream model on each gene
# By default, uses the Satterthwaite approximation for the hypothesis test

fitmm = dream( vobjDream, form, Metadata1 )

 # Examine design matrix

head(fitmm$design, 3)

# Get results of hypothesis test on coefficients of interest
#topTable( fitmm, coef='ptsd:wave', number=3 )results in 
#Error in fit$coefficients[, coef] : subscript out of bounds
#changed coef='ptsd:wave' to 'ptsd:wavew4'

topTable( fitmm, coef='ptsd:wavew4', number=3 )

#What genes are most differentially expressed?

top.table <- topTable( fitmm, coef='ptsd:wavew4', number=Inf ) 

#How many DE genes are there?

length(which(top.table$adj.P.Val < 0.05))

#Write top.table to a file

top.table$Gene <- rownames(top.table)
top.table <- top.table[,c("Gene", names(top.table)[1:6])]

write.table(top.table, file = "C:\\Users\\jdahr\\Desktop\\dream\\TopTableOutput_ptsd_wave.txt1", row.names = F, sep = "\t", quote = F)

# Compare p-values and make plot
p1 = topTable(fitDupCor, coef="ptsd", number=Inf, sort.by="none")$P.Value
p2 = topTable(fitmm, number=Inf, sort.by="none")$P.Value
P2 = topTable(fitmm, coef='ptsd', sort.by = "none", number=Inf )$p.value 

plotCompareP( p1, p2, vp$resp, dupcor$consensus)

#Gene set enrichment analysis (GSEA)

#GSEA for ptsd_cases_vs_control

df=read.csv("C:\\Users\\jdahr\\Desktop\\dream\\TopTableOutput_ptsd_cases_vs_control_final.csv")
ego=enrichGO(gene=df$ensemblid,keyType="ENSEMBL",OrgDb=org.Hs.eg.db,ont='ALL',readable=TRUE)

#GSEA for PTSD


#GSEA for PTSD_cases_vs_control_nowave

#GSEA for ptsd_wave

df1=read.csv("C:\\Users\\jdahr\\Desktop\\dream\\TopTableOutput_ptsd_wave_final.csv")
ego=enrichGO(gene=df$Gene,keyType="ENSEMBL",OrgDb=org.Hs.eg.db,ont='ALL',readable=TRUE)



#variancePartition plot
# PTSD_cases_controls
form <- ~ ptsd + resp + wave + Sex + age + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3 + cellFrac_ilr_4
vobjDream = voomWithDreamWeights( geneExpr, form, Metadata1 )
vp = fitExtractVarPartModel( vobjDream, form, Metadata1)

form <- ~ ptsd + (1|resp) + wave + Sex + age + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3 + cellFrac_ilr_4

plotVarPart( sortCols(vp))

#PTSD_wave
form <- ~ ptsd + wave + ptsd:wave + (1|resp) + Sex + age + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3 + cellFrac_ilr_4
vobjDream = voomWithDreamWeights( geneExpr, form, Metadata1 )
vp = fitExtractVarPartModel( vobjDream, form, Metadata1)
plotVarPart( sortCols(vp))


Metadata1$resp <- as.factor(Metadata1$resp)
Metadata1$wave <- as.factor(Metadata1$wave)

#trying (1 + wave | resp) sig difference between waves
form <- ~ ptsd + wave + (1 + wave | resp) + Sex + age + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3 + cellFrac_ilr_4
~ ptsd + wave + (1 + wave | resp) + Sex + age + ptsd:wave + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3 + cellFrac_ilr_4

# trying PTS symtom score
form <- ~ ptsd + (1 + wave | resp) + Sex + age + cellFrac_ilr_1 + cellFrac_ilr_2 + cellFrac_ilr_3 + cellFrac_ilr_4
