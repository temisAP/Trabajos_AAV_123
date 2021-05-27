from math import exp, log, sin, cos, tan, pi
import numpy as np
from numpy import deg2rad, rad2deg
from numpy import genfromtxt
import os

# Directorio de imágenes
figures_dir = './Figures/'
if not os.path.exists(figures_dir):
    os.makedirs(figures_dir)

# Para guardar los resultados en un csv
def save_result(name,result):
    results_dir = './Results/' # Ruta
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
    np.savetxt(results_dir+str(name)+'.csv', result, delimiter=',', fmt='%s')

# Para leer los resultados de un csv
def read_result(name):
    results_dir = './Results/' # Ruta
    if not os.path.exists(results_dir):
        print('No existe la carpeta de resultados ./Results')
    my_data = genfromtxt(results_dir+str(name)+'.csv', delimiter=',')
    return my_data

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
        self.g_fun   = self.G

        #
        self.epsilon  = 0.1
        self.verbose  = 500  #Pasos con output
        self.k_lim    = 1e5  # Máximos pasos
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
        if rfin == 6.371e6 and tfin == '':
            print('La simulación será llevada a cabo hasta llegar al nivel del mar terrestre')
            print(' Para establecer otra altidud de parada añada un argumento rfin')
            print(' >> reentrada(self,Dt,Y0,rfin = 11e6)')
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
        self.n          = -self.dy[:,0] / self.g0
        self.nmax       = self.n.max()
        print('Máxima carga alcanzada n=',self.nmax,'g')

        return self.u, self.gamma, self.r, self.theta

    # Soluciones analíticas para comparar, crea las figuras directamente
    def reentrada_analitica(self,Y0=[7e3,-10,100e3,0]):

        if Y0[2] <= self.RT and self.alt_flag == 1: Y0[2] = Y0[2] + self.RT
        z0 = ze = Y0[2]-self.RT

        rho0 = self.RHO(0)
        zs = -1/log(self.RHO(1)/rho0)

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
                z       = np.linspace(0,z0,100)
                n       = self.n

                print('Ue =',ue/1e3,'km/s')
                print('gamma_e =',rad2deg(gamma_e),'deg')
                print('z0 =',z0/1e3,'km')
                print('rho0 =',rho0,'kg/m^3')
                print('zs =',zs/1e3,'km')
            else:
                print('** Error: Corra primero la simulación numérica')
                exit()

            b = zs*rho0/(2*beta_val*sin(gamma_e))
            c = ue**2*rho0/(2*beta_val*g0)

            self.u_numerica     = u/ue
            self.u_balistica    = np.exp( b * np.exp(-z/zs) )
            self.u_planeo       = np.empty((0,1))
            for zz in z:
                self.u_planeo= np.append(self.u_planeo,(1+E_val/beta_val * self.rho_fun(zz)*RT/2)**-0.5)

            self.n_numerica     = n
            self.n_balistica    = c * np.exp(2*b*np.exp(-z/zs)-z/zs)
            self.n_planeo       = 1/E_val * (1-(self.u_planeo)**2)

            # U/Ue vs Z
            fig = plt.figure()
            plt.plot(self.u_numerica,self.altitud/1e3)
            plt.plot(self.u_balistica,z/1e3)
            plt.plot(self.u_planeo,z/1e3)
            plt.xlabel('U/Ue')
            plt.ylabel('Altitud [km]')
            plt.legend(['Numérica','Balística','Planeo'])
            plt.grid()
            plt.savefig(figures_dir+'u_analitica.pdf')
            plt.close()
            # n vs z
            fig = plt.figure()
            plt.plot(self.n_numerica,self.altitud/1e3)
            plt.plot(self.n_balistica,z/1e3)
            plt.plot(self.n_planeo,z/1e3)
            plt.xlabel('n')
            plt.ylabel('Altitud [km]')
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
        self.y  = np.array(Y)
        self.dy = np.array([0, 0, 0, 0])
        self.t  = np.array(t)

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
                Y   = Y1
                dY  = F(Y)
                t   = t+Dt
                k   = k+1
                self.y  = np.vstack([self.y,Y])
                self.dy = np.vstack([self.dy,dY])
                self.t  = np.append(self.t,t)

                Dt = Dt+0.01  # Para la siguiente vez se aumenta el paso
            elif Dt > Dt_new:             # El paso de tiempo era insuficiente y se repite la iteración
                Dt = Dt_new
            else:                         # Algún problema ha sucedido a la hora de calcular Dt_new
                Dt = 0.1

            # Criterio de parada Dt extremadamente pequeño
            if Dt < 1e-28:
                print('Imposible dar el siguiente paso')
                break

            # Criterio de parada rebote
            if Y[2] > Y0[2] + 100e3:
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

        z = r-self.RT

        g           = self.g_fun(z)
        rho         = self.rho_fun(z)
        beta_val    = self.beta(gamma)
        E_val       = self.E(gamma)

        eqns  = [-g*sin(gamma)-rho/(2*beta_val) * u**2,
        1/u * (u**2 * (rho*E_val/(2*beta_val) + cos(gamma)/r) - g*cos(gamma)),
        u*sin(gamma),
        1/r * u*cos(gamma)]

        return np.array(eqns)

    ## Funciones de gravedad y densidad

    def G(self,z):

        g0 = self.g0
        RT = self.RT
        g = g0 * (RT/(z+RT))**2

        return g

    def RHO(self,z):

        RT = self.RT
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
