# %%
%matplotlib ipympl
import numpy as np
import matplotlib.pylab as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import odeint
from scipy.optimize import *
import sys

# * En gros ce que j'ai ajouté c'est la possibilité de comparer avec une trajectoire classique (calculé avec geodesics_classical), et le calcul du redshift (qui s'affiche aussi sur le graphe)
# * J'ai aussi ajouté un calcul de l'énergie potentielle, mais je ne sais pas si c'est correct
# ! Pour l'affichage de l'énergie potentielle ce serait mieux d'avoir un truc périodique qui se répète

# Function to define the derivatives for the ODE solver (geodesic equations)
def geodesics(y, tau, L, E, M):
    t0, r, phi, r_dot = y  # phi_dot, t_dot = const
 
    dt_dtau = E/(1-2*M/r) # ! montre la "fausse" singularité en r = 2M, qui vient du problème des coordonées, mais en réalité elle est en r = 0, donc tous les calculs considéré sont pour r > 2M
    dr_dtau = r_dot
    dphi_dtau = L / r**2
    drdot_dtau = -(2* H * M /r**2) + (L**2 / r**3) - (3 * M * L**2) / r**4   #ok
    
    return [dt_dtau, dr_dtau, dphi_dtau, drdot_dtau]

def geodesics_classical(y, tau, L, E, M): # permet de comparer avec une trajectoire newtonienne classique
    t0, r, phi, r_dot = y
    
    dt_dtau = E / (1 - 2 * M / r)  # Identique
    dr_dtau = r_dot
    dphi_dtau = L / r**2
    drdot_dtau = -(2 * H * M / r**2) + (L**2 / r**3)  # Pas de terme relativiste ici
    
    return [dt_dtau, dr_dtau, dphi_dtau, drdot_dtau]

def redshift(r): # ? Est ce qu'on l'exprime comme décalage relatif ou est ce qu'on laisse les valeurs ?
    """
    Retourne une liste du redshift de la longueur d'onde d'émission au cours du temps
    """
    # G = 6.67 * 1e-11
    # c = 3 * 1e8
    
    # Calcul du redshift gravitationnel
    z = 1/np.sqrt(1 - 2*M/r) - 1
    
    return z

def compute_potential_energy(r):
    if H == 1/2:
        return - (2 * M / r - 1) * (1/2 + 1/2 * (L / r)**2)
    elif H == 0:
        #return  - (2*M/r - 1) * (1/2 * (L / r)**2)
        return L**2/(2*r**2) - M*L**2/r**3

####  Jesli od razu def w ref propre - bez conversji

# %% 
M = 3 # nombre de masse solaire M_0 = 2*1e30
H = 1/2
#dt_dtau = 1  

Lagr = H #1/2 si c'est une particule massive ou 0 si c'est un photon
L = np.sqrt(12)*M

# %% 
# * Permet de vérifier l'existence d'orbites stables
if Lagr == 1/2:
    if L >= (12)**(1/2)*M:
        print("racines existent")
        r_plus=((L**2+np.sqrt(L**4-12*(M**2)*L**2))/2/M) # rayon stable le plus éloigné de l'étoile
        r_moins=((L**2-np.sqrt(L**4-12*(M**2)*L**2))/2/M) # rayon stable le plus proche de l'étoile
        print(r_moins,r_plus)

    else:
        print("pas des racines") # cela signifie que l'étoile tombera forcément dans le trou noir
elif Lagr == 0:
    orb_stable = 3*M
# %%
# ! ##############################
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0.1, 110, 100000)  # Éviter la division par zéro
M = 3
L = 1.2 * np.sqrt(12) * M
f_x = -M/x + L**2/(2*x**2) - M*L**2/(x**3)

fig, ax = plt.subplots()
ax.plot(x, f_x)#, label=r"f(x) = -M/x + L^2/(2x^2) - M*L^2/(x^3)")

plt.title(r"Évolution de l'énergie potentielle pour $M=3$ et $L = 1.2\sqrt{12}\,M$")
plt.xlabel(r"Rayon $r$")
plt.ylabel(r"Énergie potentielle $U$")

ax.set_xlim(-1, 100)
ax.set_ylim(-0.05, 0.03)

# Supprimer les nombres des axes
ax.set_xticklabels([])
ax.set_yticklabels([])

# Garder les ticks visibles pour montrer les axes
ax.tick_params(axis='x', which='both', direction='in', length=5)
ax.tick_params(axis='y', which='both', direction='in', length=5)

ax.grid(True)
plt.show()


# ! ###############################









# %%

import numpy as np
import matplotlib.pyplot as plt
# ? euh c'est quoi ça déjà, je crois que c'était le calcul de quand l'énergie potentielle est nulle mais je sais plus pourquoi
r_pm = (-L**2/2 - np.sqrt(L**4/4 - 2*M**2*L**3)) /(-M*L) 
print(r_pm)

affich_newton = False
#69.82 million km
r = 3*M # rayon initial, peut dépendre des valeurs de rayons particuliers
# v_phi = L/r**2 # ! sert à rien pour l'instant
v_r = 0 #1*1/100 # initialisation de la vitesse radiale
t0 = 0
phi0 = 0

print(f'U = {compute_potential_energy(r)}')
print(L**2/(54*M**2))

y0 = [t0, r, phi0, v_r] 

tau = np.linspace(0, 100, 10000) # défini l'intervalle de temps choisi pour la simulation

E = (1-2*M/r)   #(1-2*M/r) A une importance seulement pour le temps propre
# ! Constant uniquement pour orbite circulaire

sol = odeint(geodesics, y0, tau, args=(L, E, M,))
sol_t0, sol_r, sol_phi, sol_r_dot = sol[:, 0], sol[:, 1], sol[:, 2], sol[:, 3]

sol_newton = odeint(geodesics_classical, y0, tau, args=(L, E, M,))
sol_t0_newton, sol_r_newton, sol_phi_newton, sol_r_dot_newton = sol_newton[:, 0], sol_newton[:, 1], sol_newton[:, 2], sol_newton[:, 3]

x_sol =  sol_r * np.cos(sol_phi)
y_sol =  sol_r * np.sin(sol_phi)

x_sol_newton =  sol_r_newton * np.cos(sol_phi_newton)
y_sol_newton =  sol_r_newton * np.sin(sol_phi_newton)

#########################################################
### Partie affichage de l'animation de la trajectoire ###
#########################################################
# %%
# Création de la figure et des axes
fig, ax1 = plt.subplots()
ax1.set_xlabel("x")
ax1.set_ylabel("y")
ax1.set_title("Évolution d'une trajectoire 2D")

# Initialisation de l'objet graphique
line, = ax1.plot([], [], lw=2, label="Trajectoire")
if affich_newton:
    line_class, = ax1.plot([], [], lw=2, label="Classique", color='orange', linestyle='--')  # Trajectoire classique
    point_class, = ax1.plot([], [], 'go', label='Position actuelle (classique)')  # Point pour classique
point, = ax1.plot([], [], 'ro')  # Point pour la position actuelle


if affich_newton:
    ax1.set_xlim(np.min(np.minimum(x_sol,x_sol_newton)), np.max(np.maximum(x_sol,x_sol_newton)))
    ax1.set_ylim(np.min(np.minimum(y_sol,y_sol_newton)), np.max(np.maximum(y_sol,y_sol_newton)))
else:
    ax1.set_xlim(np.min(x_sol), np.max(x_sol))
    ax1.set_ylim(np.min(y_sol), np.max(y_sol))
# Dessiner le trou noir comme un point noir au centre
ax1.plot(0, 0, 'o', color='black', markersize=8, label='Trou noir')  # Marqueur noir au centre
text_redshift = ax1.text(0.95, 0.05, '', transform=ax1.transAxes, fontsize=10, 
                        verticalalignment='bottom', horizontalalignment='right')

ax1.legend()

# Fonction d'initialisation
def init():
    line.set_data([], [])
    point.set_data([], [])
    text_redshift.set_text('')
    if affich_newton:
        line_class.set_data([], [])
        point_class.set_data([], [])
        return line, line_class, point, point_class, text_redshift
    return line, point, text_redshift


def update(frame):
    # Mise à jour des trajectoires
    line.set_data(x_sol[:frame], y_sol[:frame])  # Trajectoire relativiste
    point.set_data([x_sol[frame]], [y_sol[frame]])  # Point relativiste


    r_current = sol_r[frame]
    z_current = redshift(r_current)
    text_redshift.set_text(f'Redshift: {z_current:.3e}')  # Format scientifique

    if affich_newton:
        # Mise à jour des points (positions actuelles)
        line_class.set_data(x_sol_newton[:frame], y_sol_newton[:frame])    # Trajectoire classique
        point_class.set_data([x_sol_newton[frame]], [y_sol_newton[frame]])    # Point classique
        return line, line_class, point, point_class, text_redshift


    return line, point, text_redshift

ani = FuncAnimation(fig, update, frames=len(sol_t0), init_func=init, blit=True, interval=20)
#ani2 = FuncAnimation(fig, update_2, frames=len(sol_t0), init_func=init, blit=True, interval=5)
plt.show()

#E_p = (2*M/sol_r - 1)*(-1/2 - 1/2*(L/sol_r)**2)

#sys.exit()

"""
T=[]   #Partie energie a refaire
V=[]
for i in range(len(sol_r)):
    T.append(sol_r_dot[i]**2*1/2)
    #print("j",type((2*M/sol_r[i]-1)*(-1/2-1/2*L**2/sol_r[i]**2)))
    #V.append((2*M/sol_r[i]-1)*(-1/2-1/2*L**2/sol_r[i]**2)) #cours
    V.append(-M/sol_r[i]+L/(2*r**2)-M*L**2/r**3)   #exo 6 to wyglada sensownie w miarę bardziej niż z cours . DOBRE

x= np.linspace(0, r+20, 1000)
pot_r=[]

# Potentiel en fonction de r ? 
pot_r=(2*M/x-1)*(-1/2-1/2*L**2/x**2)  # STARE czy nie powinno sie anulować gdzieś? a raczej wypłaszczać, miec ekstremum?
#mozna wyznaczyc ekstremum przez taux d'accroisamnt
print("r",x)
pot_r=(M/x-L/(2*x**2)-M*L**2/x**3)
pot_r_0=(M/(x**2)-(x**2*v_phi)/(x**3)+3*M*(x**2*v_phi)**2/x**4)  #Przecoież l zależny od r oczywiscie 


#for i in x:
#print((M/(i**2)-L/(i**3)+3*M*L**2/i**4))

#plt.plot(x,pot_r,color='blue')
#plt.plot(x,pot_r_0,color='red')  #jak mozna zrobic generalna energie potencjalna jesli nie ma się L? 
plt.plot(V,color='blue')
plt.show()

EV=[]
for i in range(len(V)): 
    EV.append(T[i] - V[i]) # plus czy -
plt.plot(EV,color='blue')
#print(V)
plt.show()
#sys.exit()

plt.plot(T,color='blue')
plt.show()

# 1/2*v_r+(2*M/r-1)*(-1/2-1/2*L**2/r**2)=1/2*E**2
#plt.plot(sol_r,color='blue')
#plt.show()


### Energie cinetique a chaque instant
m=2000000

T= 1/2* m *sol[:, 1]**2 + 1/2* m*(sol[:, 1]*sol[:, 2])**2
#for i in range (0,len(T)):
   #print(sol[:, 2][i], sol[:, 3][i])
   #print(sol[:, 1][i]*sol[:, 2][i],(sol[:, 1][i]*sol[:, 2][i])**2 )  #sol 2 czeba zrobic z tego derivee
#print(sol[:, 3][-6:],sol[:, 3][:10] )
#plt.plot(T,color='blue')
#plt.show()
#print(len(T), type(T))
#sys.exit()"""

"""
plot_t="s" # s: plot d'orbite, "m" plot d'orbite et des solutions pour phi,r et vr

# Convert to Cartesian coordinates for plotting
x_sol =  sol_r * np.cos(sol_phi)
y_sol =  sol_r * np.sin(sol_phi)

if plot_t=="s":

    plt.plot(x_sol, y_sol)
    plt.xlabel('x')
    plt.ylabel('y')
    titl= f'Trajectoire de la masse. R: {str( r )}, M: {str(M )},  V_phi: {str(v_phi)}, V_r: { str( v_r)}'
    plt.scatter(x=0, y=0, c='black')
    plt.title(titl)
    plt.grid(True)
    plt.gca().set_aspect('equal', adjustable='box')


elif plot_t=="m":
    fig, axes = plt.subplots(2,2, figsize=(10,10))
    #gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1])  # More space for the top row

# First subplot - square plot
    #axes[0][0] = fig.add_subplot(gs[0, 0])

    axes[0][0].plot(x_sol, y_sol) # plot orbite
    axes[0][0].scatter(x=0, y=0, c='black')
    titl='Trajectoire de la masse. R: '+str( r )+'  M: '+ str(M )+ '  V_phi: ' +str(v_phi)+ "  V_r: "+ str( v_r)
    axes[0][0].set_title(titl)
    axes[0][0].set_xlabel("y")
    axes[0][0].set_ylabel("x")
    axes[0][0].grid(True)

    axes[0][1].plot(sol_r) 
    axes[0][1].set_title("Evolution de R au cours du temps")
    axes[0][1].ticklabel_format(axis='x', style='sci', scilimits=(4,4))
    axes[0][1].set_xlabel("temps")
    axes[0][1].set_ylabel("R")
    axes[0][1].grid(True)

    axes[1][0].plot(sol_phi) # green (g) et croix (x)
    axes[1][0].set_title("Evolution de phi au cours du temps")
    axes[1][0].ticklabel_format(axis='x', style='sci', scilimits=(4,4))
    axes[1][0].set_xlabel("temps")
    axes[1][0].set_ylabel("phi")
    axes[1][0].grid(True)

    axes[1][1].plot(sol_r_dot) # par défaut, points reliés (bleus)
    axes[1][1].set_title("Evolution du vitesse radiale au cours du temps")
    axes[1][1].ticklabel_format(axis='x', style='sci', scilimits=(4,4))
    axes[1][1].set_xlabel("temps")
    axes[1][1].set_ylabel("V_r")
    axes[1][1].grid(True)

    #donner les noms au axes!

    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.3, hspace=0.3)


elif plot_t=="en":
    pass

plt.show()
"""
######
# %%
