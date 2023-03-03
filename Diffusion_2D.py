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
import Scripts.Gammas as Gammas
import Scripts.Neighbors as Neighbors

def Cloud(p, f, v, t, triangulation = False, tt = [], implicit = False, lam = 0.5):
    """
    2D Diffusion Equation implemented in Unstructured Clouds of Points.
    
    This routine calculates an approximation to the solution of Diffusion equation in 2D using a Generalized Finite Differences scheme on unstructured clouds of points.
    
    The problem to solve is:
    
    \frac{\partial u}{\partial t}= v\nabla^2 u
    
    Input:
        p           m x 3           Array           Array with the coordinates of the nodes and the flag for boundary or inner node.
        f                           Function        Function declared with the boundary condition.
        v                           Real            Diffusion coefficient.
        t                           Integer         Number of time steps to be considered.
        triang                      Logical         Select whether or not there is a triangulation available.
                                                        True: Triangulation available.
                                                        False: No triangulation available (Default)
        tt          m x 3           Array           Array with the triangulation indexes.
        implicit                    Logical         Select whether or not use an implicit scheme.
                                                        True: Implicit scheme used.
                                                        False: Explicit scheme used (Default).
        lam                         Real            Lambda parameter for the implicit scheme.
                                                        Must be between 0 and 1 (Default: 0.5).
    
    Output:
        u_ap        m x 1           Array           Array with the approximation computed by the routine.
        u_ex        m x 1           Array           Array with the theoretical solution.
        vec         m x nvec        Array           Array with the correspondence of the 'nvec' neighbors of each node.
    """

    # Variable initialization
    m    = len(p[:,0])                                                              # The total number of nodes is calculated.
    nvec = 8                                                                        # Maximum number of neighbors for each node.
    T    = np.linspace(0,1,t)                                                       # Time discretization.
    dt   = T[1] - T[0]                                                              # dt computation.
    u_ap = np.zeros([m,t])                                                          # u_ap initialization with zeros.
    u_ex = np.zeros([m,t])                                                          # u_ex initialization with zeros.
    
    # Boundary conditions
    for k in np.arange(t):                                                          # For all time steps.
        for i in np.arange(m) :                                                     # For all the nodes.
            if p[i,2] == 1:                                                         # If the node is in the boundary.
                u_ap[i, k] = f(p[i, 0], p[i, 1], T[k], v)                           # The boundary condition is assigned.
  
    # Initial condition
    for i in np.arange(m):                                                          # For each of the nodes.
        u_ap[i, 0] = f(p[i, 0], p[i, 1], T[0], v)                                   # The initial condition is assigned.
    
    # Neighbor search for all the nodes.
    if triangulation == True:                                                       # If there are triangles available.
        vec = Neighbors.Triangulation(p, tt, nvec)                                  # Neighbor search with the proper routine.
    else:                                                                           # If there are no triangles available.
        vec = Neighbors.Cloud(p, nvec)                                              # Neighbor search with the proper routine.

    # Computation of Gamma values
    L = np.vstack([[0], [0], [2*v*dt], [0], [2*v*dt]])                              # The values of the differential operator are assigned.
    K = Gammas.Cloud(p, vec, L)                                                     # K computation with the required Gammas.
    
    # Generalized Finite Differences Method
    if implicit == False:                                                           # For the explicit scheme.
        K2 = np.identity(m) + K                                                     # Explicit formulation of K.
    else:                                                                           # For the implicit scheme.
        K2 = np.linalg.pinv(np.identity(m) - (1-lam)*K)@(np.identity(m) + lam*K)    # Implicit formulation of K.

    for k in np.arange(1,t):                                                        # For each of the time steps.
        un = K2@u_ap[:,k-1]                                                         # The new time-level is computed.
        for i in np.arange(m):                                                      # For all the nodes.
            if p[i,2] == 0:                                                         # If the node is an inner node.
                u_ap[i,k] = un[i]                                                   # Save the computed solution.
        
    # Theoretical Solution
    for k in np.arange(t):                                                          # For all the time steps.
        for i in np.arange(m):                                                      # For each of the nodes.
            u_ex[i,k] = f(p[i,0], p[i,1], T[k], v)                                  # The theoretical solution is computed.

    return u_ap, u_ex, vec

def Mesh(x, y, f, v, t, implicit = False, lam = 0.5):
    """
    2D Diffusion Equation implemented in Logically Rectangular Meshes.

    This routine calculates an approximation to the solution of Diffusion equation in 2D using a Generalized Finite Differences scheme in logically rectangular meshes.
    
    The problem to solve is:
     
    \frac{\partial u}{\partial t}= v\nabla^2 u
     
    Input:
        x           m x n           Array           Array with the coordinates in x of the nodes.
        y           m x n           Array           Array with the coordinates in y of the nodes.
        f                           Function        Function declared with the boundary condition.
        v                           Real            Diffusion coefficient.
        t                           Integer         Number of time steps considered.
        implicit                    Logical         Select whether or not use an implicit scheme.
                                                        True: Implicit scheme used.
                                                        False: Explicit scheme used (Default).
        lam                         Real            Lambda parameter for the implicit scheme.
                                                        Must be between 0 and 1 (Default: 0.5).
    
    Output:
        u_ap        m x n x t       Array           Array with the approximation computed by the routine.
        u_ex        m x n x t       Array           Array with the theoretical solution.
    """

    # Variable initialization
    m    = len(x[:,0])                                                              # The number of nodes in x.
    n    = len(x[0,:])                                                              # The number of nodes in y.
    T    = np.linspace(0,1,t)                                                       # Time discretization.
    dt   = T[1] - T[0]                                                              # dt computation.
    u_ap = np.zeros([m, n, t])                                                      # u_ap initialization with zeros.
    u_ex = np.zeros([m, n, t])                                                      # u_ex initialization with zeros.
    urr  = np.zeros([m*n, 1])                                                       # u_rr initialization with zeros.

    # Boundary conditions
    for k in np.arange(t):
        for i in np.arange(m):                                                      # For each of the nodes on the x boundaries.
            u_ap[i, 0,   k] = f(x[i, 0], y[i, 0], T[k], v)                          # The boundary condition is assigned at the first y.
            u_ap[i, n-1, k] = f(x[i, n-1], y[i, n-1], T[k], v)                      # The boundary condition is assigned at the last y.
        for j in np.arange(n):                                                      # For each of the nodes on the y boundaries.
            u_ap[0,   j, k] = f(x[0, j], y[0, j], T[k], v)                          # The boundary condition is assigned at the first x.
            u_ap[m-1, j, k] = f(x[m-1, j], y[m-1, j], T[k], v)                      # The boundary condition is assigned at the last x.
  
    # Initial condition
    for i in np.arange(m):                                                          # For each of the nodes on x.
        for j in np.arange(n):                                                      # For each of the nodes on y.
            u_ap[i, j, 0] = f(x[i, j], y[i, j], T[0], v)                            # The initial condition is assigned.

    # Computation of K with Gammas
    L  = np.vstack([[0], [0], [2*v*dt], [0], [2*v*dt]])                             # The values of the differential operator are assigned.
    K  = Gammas.Mesh(x, y, L)                                                       # K computation that include the Gammas.

    if implicit == False:                                                           # For the explicit scheme.
        K2 = np.identity(m*n) + K                                                   # Kp with an explicit formulation.
    else:                                                                           # For the implicit scheme.
        K2 = np.linalg.pinv(np.identity(m*n) \
                            - (1-lam)*K)@(np.identity(m*n) + lam*K)                 # Kp with an explicit formulation.

    # A Generalized Finite Differences Method
    for k in np.arange(1,t):                                                        # For each time step.
        for i in np.arange(m):                                                      # For each of the nodes on x.
            for j in np.arange(n):                                                  # For each of the nodes on y.
                urr[i + j*m, 0] = u_ap[i, j, k-1]                                   # urr as a row vector with all the solution.
                
        un = K2@urr                                                                 # New time level is computed.

        for i in np.arange(1,m-1):                                                  # For each of the interior nodes on x.
            for j in np.arange(1,n-1):                                              # For each of the interior nodes on y.
                u_ap[i, j, k] = un[i + j*m]                                         # u_ap values are assigned.

    # Theoretical Solution
    for k in np.arange(t):                                                          # For all the time steps.
        for i in np.arange(m):                                                      # For all the nodes on x.
            for j in np.arange(n):                                                  # For all the nodes on y.
                u_ex[i, j, k] = f(x[i, j], y[i, j], T[k], v)                        # The theoretical solution is computed.

    return u_ap, u_ex