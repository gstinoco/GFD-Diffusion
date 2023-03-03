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

def Cloud(p, vec, L):
    """
    2D Clouds of Points Gammas Computation.
     
    This function computes the Gamma values for clouds of points, and assemble the K matrix for the computations.
     
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
    K     = np.zeros([m,m])                                                         # K initialization with zeros.
    
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

def Mesh(x, y, L):
    """
    2D Logically Rectangular Meshes Gammas Computation.
     
    This function computes the Gamma values for Logically Rectangular Meshes, and assemble the K matrix for the computations.
     
    Input:
        x           m x n           Array           Array with the coordinates in x of the nodes.
        y           m x n           Array           Array with the coordinates in y of the nodes.
        L           5 x 1           Array           Array with the values of the differential operator.
     
     Output:
        K           m x m           Array           K Matrix with the computed Gammas.
    """
    # Variable initialization
    m  = len(x[:,0])                                                                # The number of nodes in x.
    n  = len(x[0,:])                                                                # The number of nodes in y.
    K  = np.zeros([(m)*(n), (m)*(n)])                                               # K initialization with zeros.

    # Gammas computation and Matrix assembly
    for i in np.arange(1,m-1):                                                      # For each of the inner nodes on x.
        for j in np.arange(1,n-1):                                                  # For each of the inner nodes on y.
            u  = np.array(x[i-1:i+2, j-1:j+2])                                      # u is formed with the x-coordinates of the stencil.
            v  = np.array(y[i-1:i+2, j-1:j+2])                                      # v is formed with the y-coordinates of the stencil
            dx = np.hstack([u[0,0] - u[1,1], u[1,0] - u[1,1], \
                            u[2,0] - u[1,1], u[0,1] - u[1,1], \
                            u[2,1] - u[1,1], u[0,2] - u[1,1], \
                            u[1,2] - u[1,1], u[2,2] - u[1,1]])                      # dx computation.
            dy = np.hstack([v[0,0] - v[1,1], v[1,0] - v[1,1], \
                            v[2,0] - v[1,1], v[0,1] - v[1,1], \
                            v[2,1] - v[1,1], v[0,2] - v[1,1], \
                            v[1,2] - v[1,1], v[2,2] - v[1,1]])                      # dy computation
            M  = np.vstack([[dx], [dy], [dx**2], [dx*dy], [dy**2]])                 # M matrix is assembled.
            M  = np.linalg.pinv(M)                                                  # The pseudoinverse of matrix M.
            YY = M@L                                                                # M*L computation.
            Gamma = np.vstack([-sum(YY), YY])                                       # Gamma values are found.
            p           = m*(j) + i                                                 # Variable to find the correct position in the Matrix.
            K[p, p]     = Gamma[0]                                                  # Gamma 0 assignation
            K[p, p-1-m] = Gamma[1]                                                  # Gamma 1 assignation
            K[p, p-m]   = Gamma[2]                                                  # Gamma 2 assignation
            K[p, p+1-m] = Gamma[3]                                                  # Gamma 3 assignation
            K[p, p-1]   = Gamma[4]                                                  # Gamma 4 assignation
            K[p, p+1]   = Gamma[5]                                                  # Gamma 5 assignation
            K[p, p-1+m] = Gamma[6]                                                  # Gamma 6 assignation
            K[p, p+m]   = Gamma[7]                                                  # Gamma 7 assignation
            K[p, p+1+m] = Gamma[8]                                                  # Gamma 8 assignation
    
    for j in np.arange(n):                                                          # For all the nodes in y.
        K[m*j, m*j] = 0                                                             # Zeros for the boundary nodes.
    
    for i in np.arange(1,m-1):                                                      # For all the nodes in x.
        p = i+(n-1)*m                                                               # Indexes for the boundary nodes.
        K[i, i] = 0                                                                 # Zeros for the boundary nodes.
        K[p, p] = 0                                                                 # Zeros for the boundary nodes.
    
    return K