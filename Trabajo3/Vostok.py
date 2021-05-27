from mision_class import mision, save_result, read_result
from numpy import rad2deg, deg2rad, cos, sin
import numpy as np

# Representaciones gráficas
import matplotlib.pyplot as plt
import os

# Directorio de imágenes
figures_dir = './Figures/Vostok/'
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)

# Reentrada de la cápsula Vostok 1
Vostok = mision()

# Coeficiente balístico y eficiencia según el ángulo de vuelo
Vostok.beta = lambda gamma: 628.0851

# Densidad ISA
from isa import isa_rho
Vostok.rho_fun = isa_rho

# Reentrada

Ue      = 7.8232e3
gamma_e = deg2rad(-5)
z0      = 315e3
theta0  = 0

Vostok.verbose = 100  # Cada cuantas iteraciones hace display
Vostok.epsilon = 0.01 # Error admisible

# Estudio de la eficiencia aerodinámica necesaria

Eficiencias = np.linspace(0,0.5,6)
#Eficiencias = [0,0.2,0.4,0.6,0.8,0.9,1]

nmax_list = []
for Eficiencia in Eficiencias:
    Vostok.E = lambda gamma: Eficiencia
    Vostok.reentrada([Ue,gamma_e,z0,theta0],rfin=Vostok.RT+7e3)
    save_result('Vostok'+str(round(Eficiencia*100))+'t',Vostok.t)
    save_result('Vostok'+str(round(Eficiencia*100))+'y',Vostok.y)
    save_result('Vostok'+str(round(Eficiencia*100))+'dy',Vostok.dy)
    nmax_list.append(Vostok.nmax)

fig = plt.figure()
legends = []
for Eficiencia in Eficiencias:
    y  = read_result('Vostok'+str(round(Eficiencia*100))+'y')
    dy = read_result('Vostok'+str(round(Eficiencia*100))+'dy')
    plt.plot(-dy[:,0]/Vostok.g0,(y[:,2]-Vostok.RT)/1e3)
    legends.append('E ='+str(round(Eficiencia,3)))
plt.xlabel('Factor de carga [g]')
plt.ylabel('Altitud [km]')
plt.grid()
plt.legend(legends)
plt.savefig(figures_dir+'Estudio_eficiencia.pdf')
plt.close()

## Estudio de las soluciones a posteriori

altura_tiempo       = 'yes'
velocidad_tiempo    = 'yes'
gamma_tiempo        = 'yes'
trayectorias        = 'yes'

if altura_tiempo == 'yes':
    fig = plt.figure()
    for Eficiencia in Eficiencias:
        t  = read_result('Vostok'+str(round(Eficiencia*100))+'t')
        y  = read_result('Vostok'+str(round(Eficiencia*100))+'y')
        dy = read_result('Vostok'+str(round(Eficiencia*100))+'dy')
        plt.plot(t,(y[:,2]-Vostok.RT)/1e3)
        legends.append('E ='+str(round(Eficiencia,3)))
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Altitud [km]')
    plt.grid()
    plt.legend(legends)
    plt.savefig(figures_dir+'Altitud.pdf')
    plt.close()

if velocidad_tiempo == 'yes':
    fig = plt.figure()
    for Eficiencia in Eficiencias:
        t  = read_result('Vostok'+str(round(Eficiencia*100))+'t')
        y  = read_result('Vostok'+str(round(Eficiencia*100))+'y')
        dy = read_result('Vostok'+str(round(Eficiencia*100))+'dy')
        plt.plot(t,(y[:,0])/1e3)
        legends.append('E ='+str(round(Eficiencia,3)))
    plt.xlabel('Tiempo [s]')
    plt.ylabel('Velocidad [km/s]')
    plt.grid()
    plt.legend(legends)
    plt.savefig(figures_dir+'Velocidad.pdf')
    plt.close()

if gamma_tiempo == 'yes':
    fig = plt.figure()
    for Eficiencia in Eficiencias:
        t  = read_result('Vostok'+str(round(Eficiencia*100))+'t')
        y  = read_result('Vostok'+str(round(Eficiencia*100))+'y')
        dy = read_result('Vostok'+str(round(Eficiencia*100))+'dy')
        plt.plot(t,rad2deg(y[:,1]))
        legends.append('E ='+str(round(Eficiencia,3)))
    plt.xlabel('Tiempo [s]')
    plt.ylabel('gamma [deg]')
    plt.grid()
    plt.legend(legends)
    plt.savefig(figures_dir+'gamma.pdf')
    plt.close()

if trayectorias == 'yes':
    fig = plt.figure()
    for Eficiencia in Eficiencias:
        t  = read_result('Vostok'+str(round(Eficiencia*100))+'t')
        y  = read_result('Vostok'+str(round(Eficiencia*100))+'y')
        dy = read_result('Vostok'+str(round(Eficiencia*100))+'dy')
        x = y[:,2]*sin(y[:,3])/1e3
        plt.plot(x,(y[:,2]-Vostok.RT)/1e3)
        legends.append('E ='+str(round(Eficiencia,3)))
    plt.xlabel('x [km]')
    plt.ylabel('Altitud [km]')
    plt.grid()
    plt.legend(legends)
    plt.savefig(figures_dir+'trayectorias.pdf')
    plt.close()


## Trayectoria de la reentrada seleccionada

trayectoria         = 'yes'

index = nmax_list.index(max([value for index,value in enumerate(nmax_list) if value < 8]))
print('')
print('***')
print('La eficiencia seleccionada es:',Eficiencias[index])
print('Esto supone una aceleración de:',nmax_list[index],'g')
print('***')
print('')

Vostok.E = lambda gamma: Eficiencias[index]
Vostok.reentrada([Ue,gamma_e,z0,theta0],rfin=Vostok.RT+7e3)

import matplotlib.image as image
import matplotlib.patches as patches

if trayectoria == 'yes':
    x = Vostok.r*sin(Vostok.theta)/1e3
    y = Vostok.r*cos(Vostok.theta)/1e3

    fig = plt.figure()
    plt.plot(x,Vostok.altitud/1e3)
    plt.grid()
    plt.xlabel('x [km]')
    plt.ylabel('Altitud [km]')
    plt.savefig(figures_dir+'Trayectoria.pdf')
    plt.close()

    fig, ax = plt.subplots()

    Tierra = plt.Circle((0, 0), Vostok.RT/1e3, color='tab:blue')
    ax.add_patch(Tierra)
    plt.plot(x,y)
    plt.grid()
    #ax.set_xlim(min(x)-1e3,max(x)+1e3)
    ax.set_ylim(min(y)-1e3,max(y)+1e3)
    ax.set_aspect('equal', 'box')
    plt.xlabel('x [km]')
    plt.ylabel('y [km]')
    plt.savefig(figures_dir+'Trayectoria_planeta.pdf')
    plt.close()

## Correlación de la reentrada seleccionada

Vostok.reentrada_analitica([Ue,gamma_e,z0,theta0])
