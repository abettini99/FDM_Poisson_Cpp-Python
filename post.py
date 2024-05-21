import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt
from src.post.readDataBinary import readDataBinary

##// ============== //##
##// Postprocessing //##
##// ============== //##
fileName = "./bin/data.bin"
x, y, u  = readDataBinary(fileName) # u[j,i], not u[i,j]

dudy,    dudx    = np.gradient(u,    y[:,0], x[0,:], edge_order=2)
du2dydy, du2dydx = np.gradient(dudy, y[:,0], x[0,:], edge_order=2)
du2dxdy, du2dxdx = np.gradient(dudx, y[:,0], x[0,:], edge_order=2)
curlgradu        = du2dydx - du2dxdy

print(curlgradu)

##// ======== //##
##// Plotting //##
##// ======== //##

# plt.contourf(x,y,u)


plt.contourf(x,y,curlgradu)
plt.show()