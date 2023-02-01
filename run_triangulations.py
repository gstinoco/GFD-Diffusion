# All the codes presented below were developed by:
#   Dr. Gerardo Tinoco Guerrero
#   Universidad Michoacana de San Nicolás de Hidalgo
#   gerardo.tinoco@umich.mx
#
# With the funding of:
#   National Council of Science and Technology, CONACyT (Consejo Nacional de Ciencia y Tecnología, CONACyT). México.
#   Coordination of Scientific Research, CIC-UMSNH (Coordinación de la Investigación Científica de la Universidad Michoacana de San Nicolás de Hidalgo, CIC-UMSNH). México
#   Aula CIMNE-Morelia. México
#
# Date:
#   January, 2023.
#
# Last Modification:
#   January, 2023.

import numpy as np
from scipy.io import loadmat
import Scripts.Errors as Errors
import Scripts.Graph as Graph
import Diffusion_2D

# Number of Time Steps
t  = 10000

# Diffusion coefficient
nu = 0.2

# Boundary conditions
# The boundary conditions are defined as
#   f = e^{-2*\pi^2vt}\cos(\pi x)cos(\pi y)

def fDIF(x, y, t, v):
    fun = np.exp(-2*np.pi**2*v*t)*np.cos(np.pi*x)*np.cos(np.pi*y)
    return fun

# Region data is loaded.
# Triangulation to work in.

regions = ['CAB','CUA','CUI','DOW','ENG','GIB','HAB','MIC','PAT','ZIR']
sizes = ['1', '2', '3']

for reg in regions:
    region = reg

    for me in sizes:
        cloud = me

        # All data is loaded from the file
        mat  = loadmat('Data/Clouds/' + region + '_' + cloud + '.mat')
        nomt = 'Results/Explicit/Triangulations/' + region + '_' + cloud + '.png'

        # Node data is saved
        p   = mat['p']
        tt  = mat['tt']
        if tt.min() == 1:
            tt -= 1

        # Poisson 2D computed in a triangulation
        u_ap, u_ex, vec = Diffusion_2D.Triangulation(p, tt, fDIF, nu, t)
        er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)
        print('The maximum mean square error in the triangulation', region, 'with size', cloud, 'is: ', er.max())
        #Graph.Error(er)
        #Graph.Cloud_Transient(p, tt, u_ap, u_ex)