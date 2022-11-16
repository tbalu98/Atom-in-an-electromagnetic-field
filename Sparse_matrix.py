import numpy as np

#This data structure records what elements of the sparse matrice is non-zero.
#For the RK4 algorithm I have to indicate that a row no longer contains non-zero elements in the matrice,
#therefore I used a 2D array to represent the indexes of the non-zero elements.
#the end of the row means that the sparse matrix row no longer contains non-zero elements
class Trip_in_sparse_matrix:
    def __init__(self, matrix):
        
        self.trip = []

        for i in range(len(matrix)):
            self.trip.append([]) 
            for j in range(len(matrix[i])): 
                if matrix[i][j]!=0:
                    self.trip[i].append((i,j))
