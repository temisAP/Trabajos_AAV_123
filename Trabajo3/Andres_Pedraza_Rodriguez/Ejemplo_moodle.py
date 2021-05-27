from mision_class import mision
from numpy import rad2deg, deg2rad, cos, sin

# Reentrada de la cápsula Viking 1
Ejemplo = mision()

# Coeficiente balístico y eficiencia según el ángulo de ataque (FLUENT)
Ejemplo.beta = lambda gamma: 64/1.605 * (-2.905e-4*abs(rad2deg(gamma))**2-0.0053*abs(rad2deg(gamma))+1.605)
Ejemplo.E    = lambda gamma: (-4.883e-4*rad2deg(gamma)**2+0.0344*rad2deg(gamma)-0.0278)/(-2.905e-4*rad2deg(gamma)**2-0.0053*rad2deg(gamma)+1.605)

# Densidad ISA
from isa import isa_rho
Ejemplo.rho_fun = isa_rho

# Reentrada

Ue      = 7.5e3
gamma_e = deg2rad(-5)
z0      = 100e3
theta0  = 0

Ejemplo.verbose = 1000  # Cada cuantas iteraciones hace display
Ejemplo.epsilon = 1     # Error admisible
Ejemplo.reentrada([Ue,gamma_e,z0,theta0])

# Representaciones gráficas
import matplotlib.pyplot as plt
import os
import numpy as np

# Directorio de imágenes
figures_dir = './Figures/Ejemplo'
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)

altura_tiempo   = 'yes'
trayectoria     = 'yes'

if altura_tiempo == 'yes':
    fig = plt.figure()
    plt.plot(Ejemplo.t,Ejemplo.altitud/1e3)
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Altitud [km]')
    plt.grid()
    plt.savefig(figures_dir+'Altitud.pdf')
    plt.close()

if trayectoria == 'yes':
    fig = plt.figure()
    plt.plot(Ejemplo.r*sin(Ejemplo.theta)/1e3,Ejemplo.r*cos(Ejemplo.theta)/1e3-Ejemplo.RT/1e3)
    plt.grid()
    plt.xlabel('x [km]')
    plt.ylabel('Altitud [km]')
    plt.savefig(figures_dir+'Trayectoria.pdf')
    plt.close()

Ejemplo.reentrada_analitica(Ejemplo.t[-1],[Ue,gamma_e,z0,theta0])
