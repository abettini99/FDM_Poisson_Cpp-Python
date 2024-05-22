import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt

import src.post.binaryData as BinaryData
import src.post.filledContour as FilledContour

##// ============== //##
##// Postprocessing //##
##// ============== //##
fileName = "./bin/data.bin"

# Get solution
x, y, u  = BinaryData.read(fileName) # u[j,i], not u[i,j]

# Get gradients of solution
dudy,    dudx    = np.gradient(u,    y[:,0], x[0,:], edge_order=2)
d2udydy, d2udydx = np.gradient(dudy, y[:,0], x[0,:], edge_order=2)
d2udxdy, d2udxdx = np.gradient(dudx, y[:,0], x[0,:], edge_order=2)

# Calculate curl(grad(u)) : should be zero everywhere (or close to it)
curlgradu = d2udydx - d2udxdy

# Calculate div(grad(u))  : should be source term everywhere (or close to it)
divgradu  = d2udxdx + d2udydy

##// ======== //##
##// Plotting //##
##// ======== //##
# solution plot
FilledContour.plot(x,y,u,         
                   levels=10,
                   title =r'$u\,\,[-]$',
                   xlabel=r'$x\,\,[-]$',
                   ylabel=r'$y\,\,[-]$')
# curl(grad(u)) plot
FilledContour.plot(x,y,curlgradu, 
                   levels=10,
                   title =r'$\nabla\times\nabla u\,\,[-]$',
                   xlabel=r'$x\,\,[-]$',
                   ylabel=r'$y\,\,[-]$')
# div(grad(u)) plot
FilledContour.plot(x,y,divgradu,  
                   levels=10,
                   title =r'$\nabla\cdot\nabla u\,\,[-]$',
                   xlabel=r'$x\,\,[-]$',
                   ylabel=r'$y\,\,[-]$')

plt.show()