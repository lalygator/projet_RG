import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint
import sys 

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


x = np.linspace(-5, 5, 100)   # Pour quoi faire?
y = np.linspace(-5, 5, 100)
x, y = np.meshgrid(x, y)


 

# Function to define the derivatives for the ODE solver (geodesic equations)
def geodesics(y, tau, L, E, M):
    t0, r, phi, r_dot= y  # phi_dot, t_dot = const
 
    dt_dtau = E*(1-2*M/r)
    dr_dtau = r_dot
    dphi_dtau = L / r**2
    drdot_dtau = -(M / r**2) + (L**2 / r**3) - (3 * M * L**2) / r**4
    
    return [dt_dtau, dr_dtau, dphi_dtau, drdot_dtau]

#ref distant

#M = 1 # correspond à une masse de 3 masse solaire M_0 = 2*1e30  # co się dzieje jesli 5
"""r0 =2 #4.815322 #0.001023
t0=0
phi0 = 0.0
v_phi = 12*1/100
v_r = 3*1/100"""



#  ref propre conversion 
"""dt_dtau=1/np.sqrt(((1-2*M/r0)-(1-2*M/r0)**(-1)*v_r**2-r0*v_phi**2))
#dt_dtau=1
r_dot0=v_r*dt_dtau
d_phi_d_tau=v_phi*dt_dtau

# quantite conserves
E=1-2*M/r0*dt_dtau
#print(((1-2*M/r0)-(1-2*M/r0)**(-1)*v_r**2-r0*v_phi**2),E,M,r0,dt_dtau)
L=r0**2*v_phi*dt_dtau"""



####  Jesli od razu def w ref propre - bez conversji
M=1
dt_dtau =1   #
r_dot0 = 0 #1*1/100
d_phi_d_tau = 5*1/100
v_phi = d_phi_d_tau
### entrer L au lieu de la vitesse? 

racine= fsolve(lambda r: M/(r**2) - (r**2*v_phi*dt_dtau)**2 / (r**3) + 3*M*(r**2*v_phi*dt_dtau)**2/  (r**4), 2)
print("a",racine)


#racine= fsolve(lambda r: M  - (v_phi*dt_dtau)**2 * r**3 + 3 * M * (v_phi*dt_dtau)**2 * r**2, 2)
#racine= fsolve(lambda r: M - 3*M*(v_phi**2)*(r**2) + 3*M*(r**2*v_phi*dt_dtau)**2, 0.1)
print("b",racine)


r=8.51580318   #22.59219582 #(7)**(1/2) # 49.09607484   #4.87223663   #22.59219582
L=r**2*v_phi*dt_dtau
print("pas de racine", L<(12)**(1/2)*M, L)
E=(1-2*M/r)*dt_dtau
print(E)
#sys.exit()
t0=0
phi0 = 0.0
r0=r
v_r=r_dot0
####



#  czy dobrze ze to v phi a nie dphi d tau
"""racine = fsolve(lambda r: -v_phi**2*r**6 - M*v_phi**2*r**5 + (6*M*v_phi**2-1-v_r**5)*r**4 + (7*M-3*M*v_r**5)*r**3 - 16*M**2*r**2 + (M/v_phi**4+12*M**3)*r - 2*M/v_phi**4, 0)
#M/(r**2) - (r**2*v_phi*(1/np.sqrt(((1-2*M/r)-(1-2*M/r)**(-1)*v_r**2-r*v_phi**2))))**2 /(r**3) + 3*M*(r**2*v_phi*(1/np.sqrt(((1-2*M/r)-(1-2*M/r)**(-1)*v_r**2-r*v_phi**2))))**2/(r**4), 2)
print(racine)
sys.exit()"""

print("param", dt_dtau, r_dot0, d_phi_d_tau, L, E )

tau = np.linspace(0, 600000, 6000000)


y0 = [t0, r0, phi0, r_dot0]


sol = odeint(geodesics, y0, tau, args=(L, E, M,))
print(sol[0:6])

# Extract the solutions
r_sol = sol[:, 1]
phi_sol = sol[:, 2]



# Convert to Cartesian coordinates for plotting
x_sol = r_sol * np.cos(phi_sol)
y_sol = r_sol * np.sin(phi_sol)

# Plot the trajectory
plt.plot(x_sol, y_sol)
plt.xlabel('x')
plt.ylabel('y')
titl='Trajectoire de la masse. R :'+str( r0 )+' M :'+ str(M )+ ' V_phi :' +str(v_phi)+ " V_r :"+ str( v_r)
plt.title(titl)
plt.grid(True)
#plt.gca().set_aspect('equal', adjustable='box')
plt.show()
plt.show()

######

### En bas: code de chatgpt peut etre faut
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# # Constants
# M = 1.0  # Mass of the black hole in natural units

# # Function to define the derivatives for the ODE solver (geodesic equations)
# def geodesics(y, tau, M):
#     r, phi, pr, pphi = y
    
#     # Derivatives
#     dr_dtau = pr
#     dphi_dtau = pphi / r**2
#     dpr_dtau = -(M / r**2) + (pphi**2 / r**3) - (3 * M * pphi**2) / r**4
#     dpphi_dtau = 0  # Angular momentum is conserved
    
#     return [dr_dtau, dphi_dtau, dpr_dtau, dpphi_dtau]

# # Initial conditions
# r0 = 10.0  # Initial radial distance
# phi0 = 0.0  # Initial angle
# pr0 = 0.0  # Initial radial momentum
# pphi0 = 4.0  # Initial angular momentum

# # Pack the initial conditions into an array
# y0 = [r0, phi0, pr0, pphi0]

# # Time array for integration (proper time steps)
# tau = np.linspace(0, 100, 1000)

# # Solve the ODEs
# sol = odeint(geodesics, y0, tau, args=(M,))

# # Extract the solutions
# r_sol = sol[:, 0]
# phi_sol = sol[:, 1]

# # Convert to Cartesian coordinates for plotting
# x_sol = r_sol * np.cos(phi_sol)
# y_sol = r_sol * np.sin(phi_sol)

# # Plot the trajectory
# plt.plot(x_sol, y_sol)
# plt.xlabel('x')
# plt.ylabel('y')
# plt.title('Particle trajectory in Schwarzschild spacetime')
# plt.grid(True)
# plt.gca().set_aspect('equal', adjustable='box')
# plt.show()
