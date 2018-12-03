rm(list=ls()) # clear workspace
# ----------------------------
# load relevant libraries
# ----------------------------
library(Matrix)
library(rmatio)
library(data.table)
# ----------------------------
# Download the expression data "URD_Dropseq_Expression_Log2TPM.txt" and the metadata "URD_Dropseq_Meta.txt" from https://portals.broadinstitute.org
# ----------------------------
#  https://portals.broadinstitute.org/single_cell/data/private/single-cell-reconstruction-of-developmental-trajectories-during-zebrafish-embryogenesis?filename=URD_Dropseq_Expression_Log2TPM.txt.gz
# ----------------------------
# For more info visit http://www.schierlab.fas.harvard.edu/urd/
# ----------------------------
# set the folder containing the dataset as current working directory
# ----------------------------
setwd("~/polybox/AAA ETH PROJECT/REAL SINGLE CELL DATA/Farrell/temp")
# ----------------------------
# Read data
# ----------------------------
full_exp_matrix<- fread("URD_Dropseq_Expression_Log2TPM.txt", header = TRUE)
gene_names <- full_exp_matrix$GENE
metadata<- fread("URD_Dropseq_Meta.txt", header = TRUE)
aa=metadata$Stage[-1]
# split by space, get the element in first index
time_info=as.numeric(sapply(strsplit(aa, "-"), "[[", 1))
exp_matrix <- as.matrix(full_exp_matrix[,-1])
#time_info <- rep(0,ncol(exp_matrix))  # There's no time info, set as 0
# ----------------------------
# Convert full matrix into a sparse matrix
# ----------------------------
sparse_exp_matrix <- as(exp_matrix[], "sparseMatrix")  

# ----------------------------
# Transpose in order to get a single_cell_dataset with rows=cells column=genes
# ----------------------------
sparse_exp_matrix<-t(sparse_exp_matrix)
cell <- sparse_exp_matrix@Dimnames[[1]]
# sparse_exp_matrix <- exp(sparse_exp_matrix)-1
# sparse_exp_matrix <- as(sparse_exp_matrix[], "sparseMatrix")  
# ----------------------------
# save dgcMatrix in file.mat
# ----------------------------
filename <- ('Farrell_expDATA.mat')

## Example how to read and write a S4 class with rmatio
## Create 'DemoS4Mat' class
setClass("DemoS4Mat",
         representation(expMatrix = "dgCMatrix",
                        cellID = "character",
                        geneNAMES = "character",
                        timeline = "numeric"
         ))
## Create a function to coerce a 'DemoS4Mat' object to a list.
setAs(from="DemoS4Mat",
      to="list",
      def=function(from)
      {
        return(list(expMatrix=from@expMatrix,
                    cellID=from@cellID,
                    geneNAMES=from@geneNAMES,
                    timeline=from@timeline
        ))
      } )
## Create a function to coerce a list to a 'DemoS4Mat' object.
setAs(from="list",
      to="DemoS4Mat",
      def=function(from)
      {
        return(new("DemoS4Mat",
                   expMatrix=from[["expMatrix"]],
                   cellID=from[["cellID"]],
                   geneNAMES=from[["geneNAMES"]],
                   timeline=from[["timeline"]]))
      } )
## Define a method to write a 'DemoS4Mat' object to a MAT file.

setMethod("write.mat",
          signature(object = "DemoS4Mat"),
          function(object,
                   filename,
                   compression,
                   version)
          {
            ## Coerce the 'DemoS4Mat' object to a list and
            ## call 'rmatio' 'write.mat' with the list.
            return(write.mat(as(object, "list"),
                             filename,
                             compression,
                             version))
          } )
## Create a new 'DemoS4Mat' object
demoS4mat <- new("DemoS4Mat",
                 expMatrix = sparse_exp_matrix,
                 cellID=cell,
                 geneNAMES = gene_names,
                 timeline = time_info
)
## Write to MAT file
write.mat(demoS4mat, filename)
## Read the MAT file
# demoS4mat.2 <- as(read.mat(filename), "DemoS4Mat")
# ## Check result
# stopifnot(identical(demoS4mat, demoS4mat.2))
# unlink(filename)
# ## End(Not run)