rm(list=ls()) # clear workspace
# ----------------------------
# load relevant libraries
# ----------------------------
library(Matrix)
library(rmatio)
library(data.table)

ngenes=50 
ncells=100

myMat<-matrix(runif(ncells*ngenes), ncol=ncells) 

# 1) Add expression data ('Large Matrix' ngenesxncells)
exp_matrix <- myMat

# 2) Add time info if available:
# numeric vector 1xncells (capture time for each cell, if not available set as a vector of zeros)
time_info=sort(sample(1:6, ncells,replace = TRUE))
#time_info <- rep(0,ncol(exp_matrix))  # There's no time info, set as 0
# ----------------------------
# Convert full matrix into a sparse matrix
# ----------------------------
sparse_exp_matrix <- as(exp_matrix[], "sparseMatrix")  

# ----------------------------
# Transpose in order to get a single_cell_dataset with rows=cells column=genes
# ----------------------------
sparse_exp_matrix<-t(sparse_exp_matrix)

# 3) Add cell ID info (optional). Character vector ncells elements
cell <- sprintf("Cell %s ", 1:ncells)

# 4) Add gene names. Character vector with ngenes elements
gene_names <- sprintf("Gene %s ", 1:ngenes)


# ----------------------------
# save dgcMatrix in file.mat
# ----------------------------
filename <- ('save_data_to_matlab.mat')


######################################################################
# Do not change here

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
