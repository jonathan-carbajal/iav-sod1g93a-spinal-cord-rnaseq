
# load libraries required for analysis
library(GEOquery)
library(limma)

#assign the data file to a data frame
dat <- read.csv("~/bioresearch-workspace/iav-sod1g93a-spinal-cord-rnaseq/RNAseq_Mouse_Gene.csv")

# pull the dataset directly from GEO
gse <- getGEO("GSE342036")
eset <- gse[[1]]

# clear strings out of the data and assign it as row names
rownames(dat) <- dat$X
dat$X <- NULL

# set the data frame as a matrix
expr <- as.matrix(dat)

# make sure that the data has no issues
dim(expr)
class(expr)
sum(is.na(expr))
sum(duplicated(rownames(expr)))

# create a metadata to describe each col name in expr
metadata <- data.frame(samples = colnames(expr), condition = c("Uninfected", "Uninfected", "Uninfected", "Infected", "Infected", "Infected"))
rownames(metadata) <- colnames(expr)

# make sure metadata aligns with expr
identical(rownames(metadata), colnames(expr))

# labels the condition and sets which one is default
metadata$condition <- factor(metadata$condition, levels = c("Uninfected", "Infected"))

# converts the conditions in metadata to numeric values and assigns it to design
design <- model.matrix(~ condition, data = metadata)

# data analysis. 
# we pass in expr to tell what to analyze and design to tell it under what
# conditions. eBayes helps to reduce the significance of outliers 
fit <- lmFit(expr, design)
fit <- eBayes(fit)

# get the results and analyze whether it is 
# statistically significant 

results <- topTable(fit, coef = "conditionInfected", sort.by = "P", number = Inf)

sum(results$adj.P.Val < 0.05)

# compare my data with the original

sum(eset$adj.P.Val < 0.05)

# visualize it

boxplot(expr)

