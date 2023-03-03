"""
All the codes presented below were developed by:
    Dr. Gerardo Tinoco Guerrero
    Universidad Michoacana de San Nicolás de Hidalgo
    gerardo.tinoco@umich.mx

With the funding of:
    National Council of Science and Technology, CONACyT (Consejo Nacional de Ciencia y Tecnología, CONACyT). México.
    Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México
    Aula CIMNE-Morelia. México

Date:
    November, 2022.

Last Modification:
    March, 2023.
"""

import numpy as np
 
def Mesh(x, y, L):
    """
    2D Mesh Gammas Computation.
    
    This function computes the Gamma values for logically rectangular meshes.
    
    Input:
        x           m x n           Array           Array with the coordinates in x of the nodes.
        y           m x n           Array           Array with the coordinates in y of the nodes.
        L           5 x 1           Array           Array with the values of the differential operator.
     
    Output:
        Gamma       m x n x 9       Array           Array with the computed gamma values.
    """

    # Variable initialization
    m     = len(x[:,0])                                                             # The number of nodes in x.
    n     = len(x[0,:])                                                             # The number of nodes in y.
    Gamma = np.zeros([m,n,9])                                                       # Gamma initialization with zeros.

    # Gammas computation
    for i in np.arange(1,m-1):                                                      # For the inner nodes in x.
        for j in np.arange(1,n-1):                                                  # For the inner nodes in y.
            dx = []                                                                 # dx initialization with zeros.
            dy = []                                                                 # dy initialization with zeros.
            for l in np.arange(i-1,i+2):                                            # For all the neighbors in x.
                for o in np.arange(j-1,j+2):                                        # For all the neighbors in y.
                    if l == i and o == j:                                           # If it is the central node.
                        pass                                                        # Do nothing.
                    else:                                                           # If it is not the central node.
                        dx.append(x[l,o] - x[i,j])                                  # Compute dx.
                        dy.append(y[l,o] - y[i,j])                                  # Compute dy.
            dx = np.array(dx)                                                       # Transform dx into an array.
            dy = np.array(dy)                                                       # Transform dy into an array.
            M = np.vstack([[dx], [dy], [dx**2], [dx*dy], [dy**2]])                  # M matrix is assembled.
            M = np.linalg.pinv(M)                                                   # The pseudoinverse of matrix M.
            YY = M@L                                                                # M*L computation.
            Gem = np.vstack([-sum(YY), YY])                                         # Gamma values are found.
            for k in np.arange(9):                                                  # For each of the Gamma values.
                Gamma[i,j,k] = Gem[k]                                               # The Gamma value is stored.
    return Gamma

def Cloud(p, vec, L):
    """
    2D Clouds of Points Gammas Computation.
     
    This function computes the Gamma values for clouds of points and triangulations and assemble the K matrix for the computations.
     
    Input:
        p           m x 3           Array           Array with the coordinates of the nodes and a flag for the boundary.
        vec         m x nvec        Array           Array with the correspondence of the 'nvec' neighbors of each node.
        L           5 x 1           Array           Array with the values of the differential operator.
     
     Output:
        K           m x m           Array           K Matrix with the computed Gammas.
    """
    # Variable initialization
    nvec  = len(vec[0,:])                                                           # The maximum number of neighbors.
    m     = len(p[:,0])                                                             # The total number of nodes.
    K     = np.zeros([m,m])
    
    # Gammas computation and Matrix assembly
    for i in np.arange(m):                                                          # For each of the nodes.
        if p[i,2] == 0:                                                             # If the node is an inner node.
            nvec = sum(vec[i,:] != -1)                                              # The total number of neighbors of the node.
            dx   = np.zeros([nvec])                                                 # dx initialization with zeros.
            dy   = np.zeros([nvec])                                                 # dy initialization with zeros.
            for j in np.arange(nvec):                                               # For each of the neighbor nodes.
                vec1  = int(vec[i, j])                                              # The neighbor index is found.
                dx[j] = p[vec1, 0] - p[i,0]                                         # dx is computed.
                dy[j] = p[vec1, 1] - p[i,1]                                         # dy is computed.
            M     = np.vstack([[dx], [dy], [dx**2], [dx*dy], [dy**2]])              # M matrix is assembled.
            M     = np.linalg.pinv(M)                                               # The pseudoinverse of matrix M.
            YY    = M@L                                                             # M*L computation.
            Gamma = np.vstack([-sum(YY), YY]).transpose()                           # Gamma values are found.
            K[i,i] = Gamma[0,0]                                                     # The corresponding Gamma for the central node.
            for j in np.arange(nvec):                                               # For each of the neighbor nodes.
                K[i, vec[i,j]] = Gamma[0,j+1]                                       # The corresponding Gamma for the neighbor node.
            
        if p[i,2] == 1:                                                             # If the node is in the boundary.
            K[i,i] = 0                                                              # Central node weight is equal to 0.
            for j in np.arange(nvec):                                               # For each of the neighbor nodes.
                K[i, vec[i,j]] = 0                                                  # Neighbor node weight is equal to 0.
    return K