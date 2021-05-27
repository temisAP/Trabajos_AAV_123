from math import exp, log, sin, cos, tan, pi
import numpy as np
from numpy import deg2rad, rad2deg
import os

# Directorio de imágenes
figures_dir = './Figures/'
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)

# Para guardar los resultados en un csv
results_dir = './Results/' # Ruta
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
def save_result(name,result):
    np.savetxt(results_dir+str(name)+'.csv', result, delimiter=',', fmt='%s')


class mision():

    # Constructor
    def __init__(self,beta= lambda gamma: 64, E = lambda gamma: 0.01):

        # Características del vehículo
        self.beta = beta    # Coeficiente balístico, por defecto es una función constante
        self.E    = E       # Eficiencia aerodinámica, por defecto es una función constante

        # Constantes

        self.RT = 6.371e6 #m        # Radio del planeta
        self.g0 = 9.81 #m/s²        # Gravedad a nivel de superficie
        self.mu = 3.986 #m³s^(-2)   # Constante gravitacional

        # Gravedad y atmósfera
        # Si no se modifican se usan las predeterminadas

        self.rho_fun = self.RHO
        self.g_fun = self.G

        #
        self.epsilon = 0.1
        self.verbose  = 500  #Pasos con output
        self.k_lim    = 1e5 # Máximos pasos
        self.alt_flag = 1
        self.plt_flag = 1

    # Ecuaciones de reentrada según esquema numérico
    def reentrada(self,Y0,rfin = 6.371e6,tfin = ''):

        print('*** Simulación en curso ***')

        # Si la distancia inicial es menor que el radio de la Tierra se asume que es un dato de altitud
        if Y0[2] <= self.RT and self.alt_flag == 1: Y0[2] = Y0[2] + self.RT
        print('')
        print('Se ha asumido que la condición inicial de distancia al centro de la Tierra es la altitud inicial')
        print(' Para desactivar esta opción modifique el booleano alt_flag')
        print(' >> reentrada.alt_flag = 0')

        # Criterio de parada
        print('')
        print('La simulación será llevada a cabo hasta llegar al nivel del mar terrestre')
        print(' Para establecer otra altidud de parada añada un argumento rfin')
        print(' >> reentrada(self,Dt,Y0,rfin = 11e3)')
        print(' Para establecer un tiempo de parada añada un argumento tfin')
        print(' >> reentrada(self,Dt,Y0,tfin = 42000)')

        print('')
        self.embedded_RK4(Y0,rfin,tfin)
        print('*** Simulación finalizada ***')

        self.u          = self.y[:,0]
        self.gamma      = self.y[:,1]
        self.r          = self.y[:,2]
        self.altitud    = self.y[:,2] - self.RT
        self.theta      = self.y[:,3]

        return self.u, self.gamma, self.r, self.theta

    # Soluciones analíticas para comparar, crea las figuras directamente
    def reentrada_analitica(self,T=1e3,Y0=[7e3,-10,100e3,0]):

        t = np.linspace(0,T,100)

        if Y0[2] <= self.RT and self.alt_flag == 1: Y0[2] = Y0[2] + self.RT
        z0 = Y0[2]-self.RT

        rho0 = self.rho_fun(0)
        zs = -1/log(self.rho_fun(1)/rho0)

        g0 = self.g0
        mu = self.mu
        RT = self.RT
        beta_val = self.beta(0)
        E_val    = self.E(0)

        # Representación gráfica
        print('')
        print('Se representará graficamente la solución numérica frente a las analítica')
        print(' Para desactivar esta opción modifique el booleano plt_flag')
        print(' >> reentrada.plt_flag = 0')
        if self.plt_flag == 1:
            import matplotlib.pyplot as plt
            if hasattr(self, 'u'):
                u       = self.u
                ue      = self.u[0]
                gamma_e = self.gamma[0]
                z       = np.linspace(0,self.altitud[0],100)
            else:
                print('** Error: Corra primero la simulación numérica')
                exit()

            self.u_numerica     = u/ue
            self.u_balistica    = np.exp( (zs*rho0)/(2*beta_val*sin(gamma_e)) * np.exp(-z/zs) )
            self.u_planeo       = ( 1-E_val/beta_val * rho0*RT/2 * np.exp(-z/zs) )**-0.5

            self.n_numerica     = np.gradient(self.u,self.t)
            c = ue**2*rho/(2*beta_val*g0)
            b = zs*rho0/(2*beta_val*sin(gamma_e))
            self.n_balistica    = c * np.exp(2*b*np.exp(-z/zs)-z0/zs)
            self.n_planeo       = 1/E_val * (1-(self.u_planeo)**2)

            # U/Ue vs Z
            fig = plt.figure()
            plt.plot(self.u_numerica,self.altitud/1e3)
            plt.plot(self.u_balistica,z/1e3)
            plt.plot(self.u_planeo,z/1e3)
            plt.xlabel('Altitud [km]')
            plt.ylabel('U/Ue')
            plt.legend(['Numérica','Balística','Planeo'])
            plt.grid()
            plt.savefig(figures_dir+'u_analitica.pdf')
            plt.close()
            # n vs z
            fig = plt.figure()
            plt.plot(self.n_numerica,self.altitud/1e3)
            plt.plot(self.n_balistica,z/1e3)
            plt.plot(self.n_planeo,z/1e3)
            plt.xlabel('Altitud [km]')
            plt.ylabel('n')
            plt.legend(['Numérica','Balística','Planeo'])
            plt.grid()
            plt.savefig(figures_dir+'n_analitica.pdf')
            plt.close()

    ## Esquemas numéricos

    def embedded_RK4(self,Y0,rfin,tfin):

        # Paso inicial
        Dt = 0.1

        # Condiciones iniciales
        Y = Y0
        t = 0

        # Ecuaciones y el solver RK4
        F = self.F
        RK4 = self.RK4

        # Vectores de estado, resultados
        self.y = np.array(Y)
        self.t = np.array(t)

        k = 1
        last_print = 0

        print('Iteracion:',0,'// Dt =',"{:.3e}".format(Dt),'// Altitud:', round((Y[2]-self.RT)/1e3,3), 'km // Tiempo:',0,'s')
        # Criterio de parada espacial
        while Y[2]>=rfin:

            # Criterio de parada temporal
            if isinstance(tfin,float):
                if t >= tfin: break

            # Cálculo del siguiente estado
            Y1 = RK4(Dt,Y)
            Y21 = RK4(Dt/2,Y)
            Y2 = RK4(Dt/2,Y21)

            # Cálculo del error
            err = Y1-Y2
            err[0] = err[0]
            err[1] = tan(err[1]/2)*1e3 # Los errores angulares se penalizan más
            err[2] = err[2]
            err[3] = tan(err[3]/2)*1e3
            mag = np.sqrt(err.dot(err))
            if mag < 1e-12:
                Dt_new = 1e12
            else:
                Dt_new = Dt*(self.epsilon/mag)**(1/(4+1))

            Dt_old = Dt

            if Dt <= Dt_new: # El paso de tiempo ha sido correcto y se usa la iteración
                Y = Y1
                t = t+Dt
                k = k+1
                self.y = np.vstack([self.y,Y])
                self.t = np.append(self.t,t)

                Dt = Dt+0.01  # Para la siguiente vez se aumenta el paso
            else:             # El paso de tiempo era insuficiente y se repite la iteración
                Dt = Dt_new

            # Criterio de parada rebote
            if Y[2] <= Y0[2] + 50e3:
                print('La cápsula ha rebotado')
                break

            # Criterio de parada por número de pasos
            if k>= self.k_lim:
                 print('Número máximo de pasos alcanzado')
                 break

            # Report de las iteraciones
            if k%self.verbose == 0 and k != last_print:
                print('Iteracion:',k,'// Dt =',"{:.3e}".format(Dt_old),'// Altitud:', round((Y[2]-self.RT)/1e3,3), 'km // Tiempo:',round(t,3),'s')
                last_print = k

    def RK4(self,Dt,Y):

        F = self.F

        # Runge Kutta 4
        K1 = F(Y)
        K2 = F(Y+Dt/2 * K1)
        K3 = F(Y+Dt/2 * K2)
        K4 = F(Y+Dt   * K3)

        return Y + Dt/6 * (K1+2*K2+2*K3+K4)


    # Ecuaciones de la dinámica de reentrada
    def F(self,Y):

        u       = Y[0]
        gamma   = Y[1]
        r       = Y[2]
        theta   = Y[3]

        g           = self.g_fun(r)
        rho         = self.rho_fun(r)
        beta_val    = self.beta(gamma)
        E_val       = self.E(gamma)

        eqns  = [-g*sin(gamma)-rho/(2*beta_val) * u**2,
        1/u * (u**2 * (rho*E_val/(2*beta_val) + cos(gamma)/r) - g*cos(gamma)),
        u*sin(gamma),
        1/r * u*cos(gamma)]

        return np.array(eqns)

    ## Funciones de gravedad y densidad

    def G(self,r):

        g0 = self.g0
        RT = self.RT
        z = r-RT
        g = g0 * (RT/(z+RT))**2

        return g

    def RHO(self,r):

        RT = self.RT
        z = r-RT
        if z <= 150e3: #m
            rho0 = 1.225 #kg/m^3
            z0 = 7.524e3 #m
        elif z> 150e3: #m
            rho0 = 3.875e-9 #kg/m^3
            z0 = 59.06e3 #m

        try:
            rho = rho0 * exp(-z/z0)
        except:
            rho = 1000 #kg/m^3
        return rho

#####################################
############## EJEMPLO ##############
#####################################

# Ejemplo básico

# viking1 = mision()
# viking1.reentrada(100,[1e3,deg2rad(-10),500e3,0])

#
#
# import matplotlib.pyplot as plt

#
# fig = plt.figure()
# plt.plot(viking1.t,viking1.altitud/1e3)
# plt.xlabel('Tiempo [s]')
# plt.ylabel('Altitud [km]')
# plt.grid()
# plt.savefig(figures_dir+'Altitud.pdf')
# plt.close()

# Modificación de coeficientes

#
# viking2 = mision(beta = lambda gamma:64-0.1*gamma, E = lambda gamma: 0.01+2*pi*gamma)
# viking2.reentrada(5,[5e3,deg2rad(0),500e3,0])

#
# import matplotlib.pyplot as plt

#
# fig = plt.figure()
# plt.plot(viking2.t,viking2.altitud/1e3)
# plt.xlabel('Tiempo [s]')
# plt.ylabel('Altitud [km]')
# plt.grid()
# plt.savefig(figures_dir+'Altitud2.pdf')
# plt.close()

# Modificación de funciones de densidad y gravedad

#
# def g(z):
    # return 9.81

#
# from isa import isa_rho as rho

#
# viking3 = mision(beta = lambda gamma:64-0.1*gamma, E = lambda gamma: 0.01+2*pi*gamma)
# viking3.rho_fun = rho
# viking3.g_fun = g

#
# viking3.reentrada(5,[5e3,deg2rad(0),500e3,0])

#
# import matplotlib.pyplot as plt

#
# fig = plt.figure()
# plt.plot(viking2.t,viking2.altitud/1e3)
# plt.xlabel('Tiempo [s]')
# plt.ylabel('Altitud [km]')
# plt.grid()
# plt.savefig(figures_dir+'Altitud3.pdf')
