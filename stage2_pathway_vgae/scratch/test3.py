import pandas as pd
import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
 
class Radbio_Dataset(torch.utils.data.Dataset):
    def __init__(self, dataset_size, train=True , ratio=0.8):
        adjacency_matrix_df = pd.read_csv('/Users/naminiyakan/Documents/BulkToST/Res-BRCA-PACSI/sim_hada.csv', header=None)

        # Convert to a PyTorch tensor
        adjacency_matrix = torch.tensor(adjacency_matrix_df.values, dtype=torch.float)

        # Step 2: Convert the adjacency matrix to edge index format
        #edge_index = adjacency_matrix.nonzero(as_tuple=False).t()
        edge_index = torch.nonzero(torch.tensor(adjacency_matrix)).t()
        
        # Step 3: Extract edge values (weights)
        #edge_values = adjacency_matrix[adjacency_matrix != 0]

        #df_edge_index = pd.DataFrame({
        #    'source_node': edge_index[0].numpy(),  # Source nodes (rows in adjacency matrix)
        #    'target_node': edge_index[1].numpy()   # Target nodes (columns in adjacency matrix)
        #})

        # DataFrame 2: Edge values (weights)
        
        #df_edge_values = pd.DataFrame({
        #    'edge_weight': edge_values.numpy()     # Weights of the edges
        #})
        
        print("Adjacency Matrix:")
        print(adjacency_matrix)
        print("\nEdge Index:")
        print(edge_index)
        #print("\nEdge Values (if weighted graph):")
        #print(df_edge_values)
        
        
        
        X = pd.read_csv('/Users/naminiyakan/Documents/BulkToST/Dataset/BRCA-PACSI/processed/Normalized_5k_BRCA_PACSI.csv',header=0,index_col=0)
        #X = X.to(torch.float32)
        #print(X.dtypes)
        
        # set training and test data size
        train_size=int(ratio*dataset_size)
        self.train=train

        self.data=(X.values.astype(np.float32),edge_index)

        if self.train:
            X=X[:train_size]
            edge_index=edge_index[:train_size]
            print(edge_index.dim())
            #edge_values=edge_values[:train_size]
            print("Training on {} examples".format(train_size))
        else:
            X=X[train_size:]
            edge_index=edge_index[train_size:]
            #edge_values=edge_values[train_size:]
            print("Testing on {} examples".format(dataset_size-train_size))
    def __getitem__(self, idx):
        "accessing one element in the dataset by index"
        edge_index = self.data[1]
        if edge_index.dim() == 3:
            print('yes')
            #edge_index = edge_index.squeeze(0)
        sample=(self.data[0][idx,...],edge_index)
        return sample
 
    def __len__(self):
        "size of the entire dataset"
        return len(self.data[0])


        
    
ratio = 0.95    
train_loader = DataLoader(dataset=Radbio_Dataset(dataset_size=2518,train=True, ratio=1), batch_size=128)
test_loader = DataLoader(dataset=Radbio_Dataset(dataset_size=2518,train=False, ratio=1),batch_size=128)
#print(train_loader)
for batch_idx in train_loader:
    batch_idx.validate()
    #print(adj_matrix[0,0:].size())
    #print(data[1][0,0:])
    #print(data[1].permute[1,0:])
    #print(data[1].size())
        
   