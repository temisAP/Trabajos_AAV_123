import matplotlib.pyplot as plt
import os
import numpy as np
from numpy import genfromtxt
from math import pi
import sympy

figures = 'no'

Nx = 100

class perfil():

    def __init__(self,S,x_lims=[],y_lims=[],r_lims=[],type=''):

        self.profile = S      #En lista para funciones a biyectivas

        if type == '2D' or (y_lims == [] and r_lims ==[]) :
            self.type_of_problem = '2D'
            self.limits = {'x':x_lims}
            self.x = x = np.linspace(self.limits['x'][0], self.limits['x'][1], Nx)

            if figures == 'yes':
                fig = plt.figure()
                ax = fig.add_subplot(111)
                ax.set_aspect('equal', 'box')

                for curve in self.profile:
                    y = curve(x)
                    ax.plot(x, y, color='tab:blue')

                ax.set_xlabel('X [m]')
                ax.set_ylabel('Y [m]').set_rotation(0)
                plt.grid()

                plt.show()
                plt.close()

        if type == 'axilsimetrico' or r_lims != []:
            self.profile = S
            self.type_of_problem = 'axilsimetrico'
            if y_lims == []: y_lims = [0,2*pi]
            self.limits = {'r':r_lims,'theta':y_lims}

            if figures == 'yes':
                fig = plt.figure()
                ax = fig.add_subplot(111, projection='3d')

                r = np.linspace(self.limits['r'][0], self.limits['r'][1], Nx)
                theta = np.linspace(self.limits['theta'][0], self.limits['theta'][1],360)

                R, P = np.meshgrid(r, theta)
                X, Y = R*np.cos(P), R*np.sin(P)

                for curve in self.profile:
                    Z = curve(R)
                    ax.plot_surface(Z, X, Y, color='tab:blue')

                ax.set_xlabel('X [m]')
                ax.set_ylabel('Y [m]')
                ax.set_zlabel('Z [m]')

                plt.show()
                plt.close()

    def normal(self,x,y=0):

        if self.type_of_problem == '2D':

            for curve in self.profile:
                xx = sympy.symbols('xx')
                u = curve(xx)
                nx = -u.diff(xx)/(-u.diff(xx)**2+1)**0.5
                ny = 1/(-u.diff(xx)**2+1)**0.5
                normal = [nx, ny]
                NX = nx.evalf(subs={xx: x})
                NY = ny.evalf(subs={xx: x})

                NORMAL = [NX,NY]

        return NORMAL

        if self.type_of_problem == 'axilsimetrico':

            for curve in self.profile:
                xx = sympy.symbols('xx')
                u = curve(xx)
                nx = -u.diff(xx)/(-u.diff(xx)**2+1)**0.5
                ny = 1/(-u.diff(xx)**2+1)**0.5
                normal = [nx, ny]
                NX = nx.evalf(subs={xx: x})
                NY = ny.evalf(subs={xx: x})

                NORMAL = [NX,NY]

        return NORMAL

def Newton(perfil):

    print(perfil.normal(1))


## Cilindro infinito (2D)
a = 'yes'
if a == 'yes':
    r = 1
    cilindro = perfil([lambda x: (r**2 - x**2)**0.5, lambda x: -(r**2 - x**2)**0.5], x_lims=[-r,r])
    Newton(cilindro)

## Esfera (3D)
b = 'no'
if b == 'yes':
    R = 1
    esfera = perfil([lambda r: (R**2-r**2)**0.5, lambda r: -(R**2-r**2)**0.5], r_lims=[0,R])


## Cofia (axilsimetrica)
c = 'no'
if c == 'yes':
    R = 2.7
    L = 10
    cofia = perfil([lambda r:(r/R)**3 * L],r_lims=[0,2.7])
