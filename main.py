import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# Generate a meshgrid for the 3D plot
x = np.linspace(-5, 5, 100)
y = np.linspace(-5, 5, 100)
x, y = np.meshgrid(x, y)
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Constants
M = 1.0  # Mass of the black hole in natural units

# Function to define the derivatives for the ODE solver (geodesic equations)
def geodesics(y, t, L, E, M):
    t, r, phi, r_dot= y
 
    dt_dtau = E*(1-2*M/r)
    dr_dtau = r_dot
    dphi_dtau = L / r**2
    drdot_dtau = -(M / r**2) + (L**2 / r**3) - (3 * M * L**2) / r**4
    
    return [dt_dtau, dr_dtau, dphi_dtau, drdot_dtau]


"""E=
L=
t0=
r0=
phi0=
r_dot0="""

y0 = [t0, r0, phi0, r_dot0]
t = np.linspace(0, tmax, N)

sol = odeint(geodesics, y0, t, args=(L,E,M,))
print(sol[0:6])

# Extract the solutions


# Convert to Cartesian coordinates for plotting
x_sol = r_sol * np.cos(phi_sol)
y_sol = r_sol * np.sin(phi_sol)

# Plot the trajectory
plt.plot(x_sol, y_sol)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Particle trajectory in Schwarzschild spacetime')
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()
plt.show()

######
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Constants
M = 1.0  # Mass of the black hole in natural units

# Function to define the derivatives for the ODE solver (geodesic equations)
def geodesics(y, tau, M):
    r, phi, pr, pphi = y
    
    # Derivatives
    dr_dtau = pr
    dphi_dtau = pphi / r**2
    dpr_dtau = -(M / r**2) + (pphi**2 / r**3) - (3 * M * pphi**2) / r**4
    dpphi_dtau = 0  # Angular momentum is conserved
    
    return [dr_dtau, dphi_dtau, dpr_dtau, dpphi_dtau]

# Initial conditions
r0 = 10.0  # Initial radial distance
phi0 = 0.0  # Initial angle
pr0 = 0.0  # Initial radial momentum
pphi0 = 4.0  # Initial angular momentum

# Pack the initial conditions into an array
y0 = [r0, phi0, pr0, pphi0]

# Time array for integration (proper time steps)
tau = np.linspace(0, 100, 1000)

# Solve the ODEs
sol = odeint(geodesics, y0, tau, args=(M,))

# Extract the solutions
r_sol = sol[:, 0]
phi_sol = sol[:, 1]

# Convert to Cartesian coordinates for plotting
x_sol = r_sol * np.cos(phi_sol)
y_sol = r_sol * np.sin(phi_sol)

# Plot the trajectory
plt.plot(x_sol, y_sol)
plt.xlabel('x')
plt.ylabel('y')
plt.title('Particle trajectory in Schwarzschild spacetime')
plt.grid(True)
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

