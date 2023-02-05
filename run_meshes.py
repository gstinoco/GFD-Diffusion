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

# Diffusion coefficient
nu = 0.2

# Boundary conditions
# The boundary conditions are defined as
#   f = e^{-2*\pi^2vt}\cos(\pi x)cos(\pi y)

def fDIF(x, y, t, v):
    fun = np.exp(-2*np.pi**2*v*t)*np.cos(np.pi*x)*np.cos(np.pi*y)
    return fun

# Names of the regions
regions = ['CAB','CUA','CUI','DOW','ENG','GIB','HAB','MIC','PAT','ZIR']

# Sizes of the clouds
sizes = ['21', '41', '81']

for reg in regions:
    regi = reg

    for me in sizes:
        mesh = me

        # Number of Time Steps
        if mesh == '21':
            t = 2000
        elif mesh == '41':
            t = 4000
        elif mesh == '81':
            t = 8000
        else:
            t = 10000

        # All data is loaded from the file
        mat = loadmat('Data/Meshes/' + regi + '_' + mesh + '.mat')
        nom = 'Results/Meshes/QME/' + regi + '_' + mesh + '.png'
        nov = 'Results/Meshes/Videos/' + regi + '_' + mesh + '.mp4'
        nop = 'Results/Meshes/Steps/' + regi + '_' + mesh + '_'

        # Node data is saved
        x  = mat['x']
        y  = mat['y']

        # Poisson 2D computed in a logically rectangular mesh
        u_ap, u_ex = Diffusion_2D.Mesh(x, y, fDIF, nu, t)
        er = Errors.Mesh_Transient(x, y, u_ap, u_ex)
        print('The maximum mean square error in the mesh', regi, 'with', mesh, 'points per side is: ', er.max())
        Graph.Error_sav(er, nom)
        Graph.Mesh_Transient_sav(x, y, u_ap, u_ex, nov)
        Graph.Mesh_Static_sav(x, y, u_ap, u_ex, nop)