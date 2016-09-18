#string_db_clusters

library(STRINGdb)
library(igraph)
library(pheatmap)
library(stringr)

#create environment
#mouse is 10090, human is 9606
string_db <- STRINGdb$new(version="10", species=9606, score_threshold=400, input_directory="")

path = "~/Google Drive/Randy = Mike/Experiments/Bionformatics -Target identification/Zinc ignorome - RE002/gene_lists/orthologs"

#make list of txt files, each txt file contains a list of genes/proteins
file.names <- dir(path, pattern = ".txt")

#iterate through each text file
for(i in 1:length(file.names)){
  mydata <- read.table(file.names[i], header=FALSE)
  
  #obtain STRING identifiers for each gene/protein
  mydata_mapped <- string_db$map(mydata, "V1", removeUnmappedRows=TRUE)
  
  #obtain protein clusters
  clustersList <- string_db$get_clusters(mydata_mapped$STRING_id)
  
  #Some genes may not fall into any clusters, so they are retained as clusters with one member. Let's get rid of those.
  number_of_clusters = 0
  for(j in 1:length(clustersList)){
    if(length(clustersList[[j]]) > 1)
      number_of_clusters = number_of_clusters + 1
  }
  #all genes of all clusters in one list
  clustersUnlist <- unlist(clustersList[1:number_of_clusters])

#obtain network of all inputted genes/proteins, create CSV of a matrix where 1 denotes a connection between two genes and 0 denotes no connection. However, the row and column names will be STRING identifiers and not genes.
  whole.network <- string_db$get_subnetwork(clustersUnlist)
  whole.network.asMatrix <- as.matrix(whole.network[])
  whole.network.asMatrix <- whole.network.asMatrix[match(clustersUnlist, rownames(whole.network.asMatrix)), ]
  whole.network.asMatrix <- whole.network.asMatrix[, match(clustersUnlist, colnames(whole.network.asMatrix))]
  write.csv(whole.network.asMatrix, file = paste(substr(file.names[i],1,5),"_string_network.csv", payload_id=NULL))

#obtain gene names from STRING identifiers. These you can copy paste into the CSV to replace the STRING ids
  string_codes_network <- colnames(whole.network.asMatrix)
  string_codes_data <- unlist(mydata_mapped['STRING_id'])
  genes_data <- unlist(mydata_mapped['V1'])
  indices <- c()

  for(code in string_codes_network){
    idx = match(code, string_codes_data)
    indices <- c(indices, idx)
  }

  gene_list <- c()
  for(idx in indices){
    gene_list <- c(gene_list, genes_data[idx])
  }

  write(gene_list, file = paste(substr(file.names[i],1,5),"_string_network_gene_list.txt", payload_id=NULL))

  #Export list of cluster lengths for use in Python "number_of_connections" file
  list_of_lengths <- c()
  for(k in seq(1:number_of_clusters)){
    length_of_cluster <- length(clustersList[[k]])
    list_of_lengths <- c(list_of_lengths, length_of_cluster)
  }
  list_of_lengths <- as.character(list_of_lengths)
  index_of_third_underscore_of_filename = which(strsplit(file.names[i], "")[[1]]=="_")
  file_to_write_to <- file(paste(substr(file.names[i],1,index_of_third_underscore_of_filename[3]-1),"cluster_lengths.txt"))
  writeLines(list_of_lengths, file_to_write_to, sep = "\n")
  
  
  #Export images of clusters
  for(k in seq(1:number_of_clusters)){
    string_db$get_png(clustersList[[k]], required_score=NULL, network_flavor="evidence",
                      file=paste('plot_',k,substr(file.names[i],1,5),'.png'), payload_id=NULL)
    }
}