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
v = 0.2

# Names of the regions
regions = ['CAB','CUA','CUI','DOW','ENG','GIB','HAB','MIC','PAT','ZIR']

# Sizes of the clouds
sizes = ['1', '2', '3']

for reg in regions:
    regi = reg

    for me in sizes:
        cloud = me

        # Number of Time Steps
        if cloud == '1':
            t = 2000
        elif cloud == '2':
            t = 4000
        elif cloud == '3':
            t = 8000
        else:
            t = 10000

        # Boundary conditions
        # The boundary conditions are defined as
        #   f = e^{-2*\pi^2vt}\cos(\pi x)cos(\pi y)

        def fDIF(x, y, t, v):
            fun = np.exp(-2*np.pi**2*v*t)*np.cos(np.pi*x)*np.cos(np.pi*y)
            return fun

        # All data is loaded from the file
        mat = loadmat('Data/Clouds/' + regi + '_' + cloud + '.mat')
        nom = 'Results/Clouds/QME/' + regi + '_' + cloud + '.png'
        nov = 'Results/Clouds/Videos/' + regi + '_' + cloud + '.mp4'
        nop = 'Results/Clouds/Steps/' + regi + '_' + cloud + '_'

        # Node data is saved
        p   = mat['p']
        tt  = mat['tt']
        if tt.min() == 1:
            tt -= 1

        # Poisson 2D computed in an unstructured cloud of points
        u_ap, u_ex, vec = Diffusion_2D.Cloud(p, fDIF, v, t)
        er = Errors.Cloud_Transient(p, vec, u_ap, u_ex)
        print('The maximum mean square error in the unstructured cloud of points', regi, 'with size', cloud, 'is: ', er.max())
        Graph.Error_sav(er,nom)
        Graph.Cloud_Transient_sav(p, tt, u_ap, u_ex, nov)
        Graph.Cloud_Static_sav(p, tt, u_ap, u_ex, nop)