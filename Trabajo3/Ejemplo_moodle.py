from mision_class import mision
from numpy import rad2deg, deg2rad, cos, sin

# Reentrada de la cápsula Viking 1
Ejemplo = mision()

# Coeficiente balístico y eficiencia según el ángulo de vuelo
Ejemplo.beta = lambda gamma: 500
Ejemplo.E    = lambda gamma: 0.25

# Densidad ISA
from isa import isa_rho
Ejemplo.rho_fun = isa_rho

# Reentrada

Ue      = 7.5e3
gamma_e = deg2rad(-5)
z0      = 100e3
theta0  = 0

Ejemplo.verbose = 1000  # Cada cuantas iteraciones hace display
Ejemplo.epsilon = 0.01  # Error admisible
Ejemplo.reentrada([Ue,gamma_e,z0,theta0])

# Representaciones gráficas
import matplotlib.pyplot as plt
import os
import numpy as np

# Directorio de imágenes
figures_dir = './Figures/Ejemplo/'
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)

altura_tiempo       = 'yes'
velocidad_tiempo    = 'yes'
gamma_tiempo        = 'yes'
trayectoria         = 'yes'

if altura_tiempo == 'yes':
    fig = plt.figure()
    plt.plot(Ejemplo.t,Ejemplo.altitud/1e3)
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Altitud [km]')
    plt.grid()
    plt.savefig(figures_dir+'Altitud.pdf')
    plt.close()

if velocidad_tiempo == 'yes':
    fig = plt.figure()
    plt.plot(Ejemplo.t,Ejemplo.u/1e3)
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Velocidad [km/s]')
    plt.grid()
    plt.savefig(figures_dir+'Velocidad.pdf')
    plt.close()

if gamma_tiempo == 'yes':
    fig = plt.figure()
    plt.plot(Ejemplo.t,rad2deg(Ejemplo.gamma))
    plt.xlabel('Tiempo [s]')
    plt.ylabel('gamma [deg]')
    plt.grid()
    plt.savefig(figures_dir+'gamma.pdf')
    plt.close()

if trayectoria == 'yes':
    x = Ejemplo.r*sin(Ejemplo.theta)/1e3
    y = Ejemplo.r*cos(Ejemplo.theta)/1e3

    fig = plt.figure()
    plt.plot(x,y-Ejemplo.RT/1e3)
    plt.grid()
    plt.xlabel('x [km]')
    plt.ylabel('Altitud [km]')
    plt.savefig(figures_dir+'Trayectoria.pdf')
    plt.close()

    fig, ax = plt.subplots()
    Tierra = plt.Circle((0, 0), Ejemplo.RT/1e3, color='b')
    ax.add_patch(Tierra)
    plt.plot(x,y)
    plt.grid()
    ax.set_aspect('equal', 'box')
    ax.set_xlim(min(x)-1e3,max(x)+1e3)
    ax.set_ylim(min(y)-1e3,max(y)+1e3)
    plt.xlabel('x [km]')
    plt.ylabel('y [km]')
    plt.savefig(figures_dir+'Trayectoria_planeta.pdf')
    plt.close()

Ejemplo.reentrada_analitica([Ue,gamma_e,z0,theta0])
