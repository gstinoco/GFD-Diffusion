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
from scipy.io import loadmat
import Scripts.Errors as Errors
import Scripts.Graph as Graph
import Diffusion_2D

# Diffusion coefficient
v = 0.2

# Names of the regions
regions = ['CAB','CUA','CUI','DOW','ENG','GIB','HAB','MIC','PAT','ZIR']

# Sizes of the Triangulations
sizes = ['1', '2', '3']

# Boundary conditions
# The boundary conditions are defined as
#   f = e^{-2*\pi^2vt}\cos(\pi x)cos(\pi y)

def fDIF(x, y, t, v):
    fun = np.exp(-2*np.pi**2*v*t)*np.cos(np.pi*x)*np.cos(np.pi*y)
    return fun

for reg in regions:
    regi = reg

    for me in sizes:
        cloud = me

        # Number of Time Steps
        if cloud == '1':
            t = 1000
        elif cloud == '2':
            t = 4000
        elif cloud == '3':
            t = 16000
        else:
            t = 32000

        # All data is loaded from the file
        mat = loadmat('Data/Clouds/' + regi + '_' + cloud + '.mat')

        # Node data is saved
        p   = mat['p']
        tt  = mat['tt']
        if tt.min() == 1:
            tt -= 1

        # Poisson 2D computed in an unstructured cloud of points
        u_ap, u_ex, vec = Diffusion_2D.Cloud(p, fDIF, v, t, implicit = False, triangulation = True, tt = tt, lam = 0.5)

        # Error computation
        er = Errors.Cloud(p, vec, u_ap, u_ex)
        print(regi, 'size', cloud, '. Explicit scheme: ', er.max())

        # Results
        nom = 'Results/Triangulations/Explicit/QME/' + regi + '_' + cloud + '.png'
        nov = 'Results/Triangulations/Explicit/Videos/' + regi + '_' + cloud + '.mp4'
        nop = 'Results/Triangulations/Explicit/Steps/' + regi + '_' + cloud + '_'
        Graph.Error_sav(er,nom)
        Graph.Cloud_Transient_sav(p, tt, u_ap, u_ex, nov)
        Graph.Cloud_Static_sav(p, tt, u_ap, u_ex, nop)
    
for reg in regions:
    regi = reg

    for me in sizes:
        cloud = me

        # Number of Time Steps
        if cloud == '1':
            t = 1000
        elif cloud == '2':
            t = 4000
        elif cloud == '3':
            t = 16000
        else:
            t = 32000

        # All data is loaded from the file
        mat = loadmat('Data/Clouds/' + regi + '_' + cloud + '.mat')

        # Node data is saved
        p   = mat['p']
        tt  = mat['tt']
        if tt.min() == 1:
            tt -= 1

        # Poisson 2D computed in an unstructured cloud of points
        u_ap, u_ex, vec = Diffusion_2D.Cloud(p, fDIF, v, t, implicit = True, triangulation = True, tt = tt, lam = 0.5)

        # Error computation
        er = Errors.Cloud(p, vec, u_ap, u_ex)
        print(regi, 'size', cloud, '. Implicit scheme: ', er.max())

        # Results
        nom = 'Results/Triangulations/Implicit/QME/' + regi + '_' + cloud + '.png'
        nov = 'Results/Triangulations/Implicit/Videos/' + regi + '_' + cloud + '.mp4'
        nop = 'Results/Triangulations/Implicit/Steps/' + regi + '_' + cloud + '_'
        Graph.Error_sav(er,nom)
        Graph.Cloud_Transient_sav(p, tt, u_ap, u_ex, nov)
        Graph.Cloud_Static_sav(p, tt, u_ap, u_ex, nop)