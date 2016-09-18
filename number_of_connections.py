"""
Created on Tue Jul 12 13:15:48 2016
@author: Randall Ellis
"""

#The purpose of this script is to import spreadsheets of the connections between all genes in a protein interaction network along with the cluster lengths of the clusters in that network, and sorting the spreadsheet into clusters to determine how many different clusters each gene connects to. In these initial spreadsheets, all genes in the network make up the first column and the first row. A matrix is displayed of 0s and 1s denoting whether any pair of genes is connected or not.

def number_of_connections(wdir):
    import pandas as pd
    import numpy as np
    import os, os.path

    path = os.listdir(wdir)
    network_files = []
    network_length_files = []
    
    #make lists of file names for cluster lengths (txt) and network files (csv)
    for f,i in zip(path,range(len(path))):
        if f[-3:] == 'txt':
            network_length_files.append(f)
        elif f[-3:] == 'csv':
            network_files.append(f)

#Iterate through network files and cluster length files
    for network_file, network_length_file in zip(range(len(network_files)), range(len(network_length_files))):
        #Import network file as pandas dataframe, import cluster lengths as a list of integers
            network_df = pd.read_csv(network_files[network_file], skipinitialspace=True, index_col=0)
            network_lengths = [length.rstrip('\n') for length in open(network_length_files[network_length_file])]
            network_lengths = [int(length) for length in network_lengths]
            
            #Initialize cluster connections matrix with all zero values
            network_matrix = np.zeros((len(network_df), len(network_lengths)))
            
            #Iterate through each gene in the network
            for gene in range(len(network_df)):
                row = network_df.ix[gene]
                start_index = 0
                end_index = int(network_lengths[0])
                #by using the cluster lengths, sums are taken for the connections from each gene to all genes in a single cluster
                for network_length in range(len(network_lengths)):
                    network_matrix[gene, network_length] = sum(row[start_index:end_index])
                    start_index=int(end_index)
                    if end_index < sum(network_lengths):
                        end_index+=network_lengths[network_length+1]
                    else:
                        break

#Make clusters column header
            clusters = []
            for i in range(len(network_lengths)):
                clusters.append(i)

#Convert numpy matrix into Pandas dataframe
            pandas_matrix = pd.DataFrame(data=network_matrix[0:,0:],
                                            index=network_df.index.values,
                                            columns=clusters)

#At this point, your matrix/Dataframe shows integers for all connections between each gene and each cluster. To easily get the number of unique clusters each gene connects to, let's convert all non-zero values to ones.
            for i in np.nditer(pandas_matrix, op_flags=['readwrite']):
                if i > 0:
                    i[...] = 1
            sums = []
            for i in range(len(pandas_matrix)):
                matrix_row = pandas_matrix.ix[i]
                sums.append(sum(matrix_row[:]))
            
#write your file to CSV
            pandas_matrix.to_csv(network_files[network_file][0:5] + "_cluster_connections_matrix.csv", delimiter=",")
