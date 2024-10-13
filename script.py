import numpy as np
import math
import matplotlib.pyplot as plt
import numba
from numba import jit
# from IPython.display import HTML
# from matplotlib.ticker import MultipleLocator
from time import time
import datetime
import os

'''

En este script hay dos simulaciones, la primera es para la función de onda sin campo magnético y la segunda es para la función de onda con campo magnético.
Se dan Nt = 1e6 pasos temporales para simular 0.1s. He calculado que guardando un frame por cada 1000 pasos podemos hacer un video de 20 segundos de 50fps,
que se supone que es buena calidad. 

'''

print("Comienza la ejecución del script...")

# Para guardar los 1000 frames en una carpeta
folder_name = "solutions"
os.makedirs(folder_name, exist_ok=True)


# physical quantities
e_charge = -1.602176634e-19
hbar = 1.05457182e-30
e_mass = 9.1093837e-31

# lattice size in space and time
L = 6 # in centimeter
T = 0.1 # in seconds
x_0 = 3 # in centimeter
y_0 = 3 # in centimeter

# lattice quantities
Nx = 400  # Number of x-axis points
Ny = 400 # Number of y-axis points
Nt = 1000000 # Number of temporal-axis points

# lattice spacing
a_t = T / Nt
a_x = L / Nx
a_y = L / Ny

# lattice initialization
vec_x = np.linspace(-x_0, L - x_0, Nx)
vec_y = np.linspace(-y_0, L - y_0, Ny)

x, y = np.meshgrid(vec_x, vec_y)

# Gaussian wave packet initialization
n_kx = 20 # x momentum
n_ky = 0 # y momentum

# Initial wave packet on lattice
psi_0 = np.zeros([2, Ny, Nx],dtype=np.complex64) # first row of two integer is used for time evolution
psi_0[0] = np.exp(-((y) ** 2 / (2 * 1 ** 2) + ( x + 1) ** 2 / (2 * 0.5 ** 2)) + 2 * math.pi / L * 1j * n_kx * x)

# Normalization to 1
normal = np.sum(np.absolute(psi_0[0]) ** 2) * a_x * a_y 
psi_0[0] = psi_0[0] / np.sqrt(normal)

# psi_0 plot
plt.imshow((np.absolute(psi_0[0]) ** 2), cmap='inferno', origin='lower', extent=[- x_0, L - x_0,- y_0, L - y_0])
plt.colorbar(label=f'Probability density')
plt.xlabel('x [cm]')
plt.ylabel('y [cm]')
plt.title('Initial wave packet')

plt.savefig("initial_wave_packet.png")
plt.close()


print("Se ha guardado la imagen del paquete de ondas inical en 'initial_wave_packet.png'")


# Potential function 
V = np.zeros([Ny,Nx])

# Positioning of the wall on x axis
N_min_wall = 190
N_max_wall = 210
# to cm
x_min_wall = N_min_wall * a_x + x_0
x_max_wall = N_max_wall * a_x + x_0
wall_thick = round(np.absolute(x_max_wall-x_min_wall),2)

# Positioning of slits (1 and 2) width on y axis
N_left_s1 = 160
N_right_s1 = 180
N_left_s2 = 220
N_right_s2 = 240

left_s1 = N_left_s1 * a_y - y_0
right_s1 = N_right_s1 * a_y - y_0
left_s2 = N_left_s2 * a_y - y_0
right_s2 = N_right_s2 * a_y - y_0


# Slits width (l) and distance between their center (d)
l = round(np.absolute(left_s1-right_s1),2)
d = round(np.absolute((left_s1 + right_s1)/2 - (left_s2 + right_s2)/2 ),2)


# Assigning a value to the potential where the wall is placed
for i in range(N_min_wall, N_max_wall):
    for j in range(0, N_left_s1):
        V[j, i] = 1
    for j in range(N_right_s1, N_left_s2):
        V[j, i] = 1
    for j in range(N_right_s2, Ny):
        V[j, i] = 1

# Check slits and wall position
print("Slits lenght:" , l)
print("Slits distance:" , d)
print("Wall thickness:" , wall_thick)

# l = 0.3 cm y d = 0.9 cm, son los valores con los que en el proyecto se obtienen los mejores resultados!

#-----------------------------------------------------

# Magnetic vector potential

# Position of the singularity
N_centre_x = 200
centre_x = N_centre_x * a_x - x_0

# Vector potential function
def vector_pot(x, y):
    r_squared = (x-centre_x)**2 + (y)**2
    mask = (r_squared < 0.01)
    
    u = np.where(mask, 0, -(y) / r_squared)
    v = np.where(mask, 0, (x-centre_x) / r_squared)
    return u, v

# Vector potential component on lattice
Ax, Ay = vector_pot(x, y)

# Physical quantites
B = 8e-8 # Magnetic field
R = 0.01 # Solenoid radius
phi = math.pi * R ** 2 * abs(B) # Magnetic flux

# Coefficient of the finite difference formula (see .pdf)
gamma1 = hbar / e_mass
gamma2 = phi * e_charge / (2 * math.pi * e_mass)
gamma3 = (phi * e_charge) ** 2 / (8 * math.pi ** 2 * e_mass * hbar)

C1_x = gamma1 * a_t / a_x ** 2
C1_y = gamma1 * a_t / a_y ** 2
C2_x = gamma2 * a_t / a_x
C2_y = gamma2 * a_t / a_y 
C3 = gamma3 * a_t

# Check
print("Numerical coefficients:" ,C1_x,C2_x,C3)

# Se espera que sean: Numerical coefficients: 0.005145228284397183 -2.9313666851029666e-05 8.35036869849819e-08.

# potential plot
plt.imshow(V+Ax, cmap='inferno', origin='lower', extent=[- x_0, L - x_0,- y_0, L - y_0])
plt.colorbar(label=f'Value of the potential')
plt.xlabel('x [cm]')
plt.ylabel('y [cm]')
plt.title('Potential function')
plt.savefig("potential.png")
plt.close()

print("Se ha guardado la imagen de la función potencial en 'potential.png'")

# Finite difference formula
'''
# Esta es la función original. Como nos interesa escribir los resultados en un archivo de texto, 
la he modificado para que devuelva el resultado en vez de guardarlo en una variable. "np.save" no se puede
utilizar dentro de una función optimizada con numba... Así que lo haremos en dos pasos. 

@numba.jit("c16[:,:,:](c16[:,:,:],float32,float32,float32,float32,float32)", nopython=True, nogil=True)
def compute_psi(psi,C1_x,C1_y,C2_x,C2_y,C3):
    for t in range(0, Nt-1):
        for i in range(0, Ny):
            ip=(Ny+i+1)%(Ny)
            im= (Ny+i-1)%(Ny)
            for k in range(0, Nx):
                kp=(Nx+k+1)%(Nx)
                km= (Nx+k-1)%(Nx)
                if V[i][k] == 0: # solution is computed outside the wall   
                    psi[1][i][k] = psi[0][i][k] * (1 - 1j * a_t * V[i][k] - 1j * C1_x - 1j * C1_y - 1j * C3 * (Ax[i][k] ** 2 + Ay[i][k] ** 2) ) \
                 + psi[0][ip][k] * (1j * C1_y / 2. + C2_y / 2. * Ay[i][k]) \
                 + psi[0][im][k] * (1j * C1_y / 2. - C2_y / 2. * Ay[i][k]) \
                 + psi[0][i][kp] * (1j * C1_x / 2. + C2_x / 2. * Ax[i][k]) \
                 + psi[0][i][km] * (1j * C1_x / 2. - C2_x / 2. * Ax[i][k])
                else: # we are inside the wall, solution is set to zero
                     psi[1][i][k] = 0
                        
        psi[0] = psi[1]

    return psi
'''

@numba.jit("c16[:,:,:](c16[:,:,:],float32,float32,float32,float32,float32)", nopython=True, nogil=True)
def compute_psi(psi,C1_x,C1_y,C2_x,C2_y,C3):
    for i in range(0, Ny):
        ip=(Ny+i+1)%(Ny)
        im= (Ny+i-1)%(Ny)
        for k in range(0, Nx):
            kp=(Nx+k+1)%(Nx)
            km= (Nx+k-1)%(Nx)
            if V[i][k] == 0: # solution is computed outside the wall   
                psi[1][i][k] = psi[0][i][k] * (1 - 1j * a_t * V[i][k] - 1j * C1_x - 1j * C1_y - 1j * C3 * (Ax[i][k] ** 2 + Ay[i][k] ** 2) ) \
             + psi[0][ip][k] * (1j * C1_y / 2. + C2_y / 2. * Ay[i][k]) \
             + psi[0][im][k] * (1j * C1_y / 2. - C2_y / 2. * Ay[i][k]) \
             + psi[0][i][kp] * (1j * C1_x / 2. + C2_x / 2. * Ax[i][k]) \
             + psi[0][i][km] * (1j * C1_x / 2. - C2_x / 2. * Ax[i][k])
            else: # we are inside the wall, solution is set to zero
                 psi[1][i][k] = 0
                        
    psi[0] = psi[1]

    return psi



print("Ahora empiezan las simulaciones largas ---------------------------------------")

# Comienzo
t_start = time()
print(f"{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')} - Empieza la simulación con campo magnético")

# solution is computed with magnetic field
# Aquí guardamos los frames que nos interesan!
psi = psi_0.astype(complex)
k = 0 # Para numerar los archivos.
for t in range(0, Nt-1):
    psi = compute_psi(psi,C1_x,C1_y,C2_x,C2_y,C3)
    if t % 1000 == 0 or t>=Nt-2:
        np.save(f'solutions/psi_B_05mu_400pt_t={k}.npy', np.abs(psi[1])**2)
        k += 1

# solution_shift = np.abs(compute_psi(psi_0.astype(complex),C1_x,C1_y,C2_x,C2_y,C3)[1])**2
solution_shift = np.abs(psi[1])**2

# Fin
t_end = time()

print(f"{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')} - Fin de la simulación con campo magnético")
print("Tiempo transcurrido (min): ", (t_start-t_end)/60)

print("--------------------------------------------------------------------------------------------")

# Comienzo
t_start = time()

print(f"{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')} - Empieza la simulación sin campo magnético")

# solution is computed with no magnetic field
# Aquí guardamos los frames que nos interesan!
psi = psi_0.astype(complex)
k = 0 # Para numerar los archivos.
for t in range(0, Nt-1):
    psi = compute_psi(psi,C1_x,C1_y,0,0,0)
    if t % 1000 == 0 or t>=Nt-2:
        np.save(f'solutions/psi_B_0_400pt_t={k}.npy', np.abs(psi[1])**2)
        k += 1

# solution = np.abs(compute_psi(psi_0.astype(complex),C1_x,C1_y,0,0,0)[1])**2
solution = np.abs(psi[1])**2

# Fin
t_end = time()

print(f"{datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')} - Fin de la simulación sin campo magnético")
print("Tiempo transcurrido (min): ", (t_start-t_end)/60)

print("Final de las simulaciones largas! --------------------------------------------------------")

# Plot of the solutions
fig, axs = plt.subplots(1, 2, figsize=(12, 6)) 

# first subplot
axs[0].imshow((solution), cmap='inferno', origin='lower', extent=[-3, 3, -3, 3])
axs[0].set_title(r'Solution with $B=0$')
axs[0].set_xlabel('x [cm]')
axs[0].set_ylabel('y [cm]')
axs[0].set_aspect('equal')
axs[0].grid(True)
axs[0].set_xlim([-3, 3])
axs[0].set_ylim([-3, 3])
axs[0].set_xticks([-3, -2, -1, 0, 1, 2, 3])

# Second subplot
axs[1].imshow((solution_shift), cmap='inferno', origin='lower', extent=[-3, 3, -3, 3])
axs[1].set_title(r'Solution with $B\neq0$')
axs[1].set_xlabel('x [cm]')
axs[1].set_ylabel('y [cm]')
axs[1].set_aspect('equal')
axs[1].grid(True)
axs[1].set_xlim([-3, 3])
axs[1].set_ylim([-3, 3])
axs[1].set_xticks([-3, -2, -1, 0, 1, 2, 3])

plt.tight_layout()
plt.savefig("Solutions_heatmap.png")

print("Se ha guardado la imagen de las soluciones en 'Solutions_heatmap.png'")

# print the solutions
with open("solution_B_0_400pt.txt", "w") as file:
    for row in solution:
        file.write(f"[{', '.join(map(str, row))}]\n")
        
with open("solution_B_05mu_400pt.txt", "w") as file:
    for row in solution_shift:
        file.write(f"[{', '.join(map(str, row))}]\n")

print("Se han guardado los resultados en 'solution_B_0_400pt.txt' y 'solution_B_05mu_400pt.txt'")