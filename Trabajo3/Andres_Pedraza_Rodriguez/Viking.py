from mision_class import mision
from numpy import rad2deg, deg2rad, cos, sin


# Reentrada de la cápsula Viking 1
Viking1 = mision()

# Coeficiente balístico y eficiencia según el ángulo de ataque (FLUENT)
Viking1.beta = lambda gamma: 64/1.605 * (-2.905e-4*abs(rad2deg(gamma))**2-0.0053*abs(rad2deg(gamma))+1.605)
Viking1.E    = lambda gamma: (-4.883e-4*rad2deg(gamma)**2+0.0344*rad2deg(gamma)-0.0278)/(-2.905e-4*rad2deg(gamma)**2-0.0053*rad2deg(gamma)+1.605)

# Densidad y gravedad marcianas
from mars_atmosphere import mars_atm_rho as mars_rho
Viking1.rho = mars_rho
Viking1.g0  = 3.721     #m/s²
Viking1.RT  = 3.3895e6  #m
Viking1.mu  = 4.282e13  #m³s^(-2)

# Reentrada

Ue      = 0.5e3
gamma_e = deg2rad(-15)
z0      = 300e3
theta0  = 0

Viking1.verbose = 1000  # Cada cuantas iteraciones hace display
Viking1.epsilon = 1     # Error admisible
Viking1.reentrada([Ue,gamma_e,z0,theta0],rfin=Viking1.RT)

# Representaciones gráficas
import matplotlib.pyplot as plt
import os
import numpy as np

# Directorio de imágenes
figures_dir = './Figures/Viking1'
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)

altura_tiempo   = 'yes'
trayectoria     = 'yes'

if altura_tiempo == 'yes':
    fig = plt.figure()
    plt.plot(Viking1.t,Viking1.altitud/1e3)
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Altitud [km]')
    plt.grid()
    plt.savefig(figures_dir+'Altitud.pdf')
    plt.close()

if trayectoria == 'yes':
    fig = plt.figure()
    plt.plot(Viking1.r*sin(Viking1.theta)/1e3,Viking1.r*cos(Viking1.theta)/1e3-Viking1.RT/1e3)
    plt.grid()
    plt.xlabel('x [km]')
    plt.ylabel('Altitud [km]')
    plt.savefig(figures_dir+'Trayectoria.pdf')
    plt.close()

Viking1.reentrada_analitica(Viking1.t[-1],[Ue,gamma_e,z0,theta0])
