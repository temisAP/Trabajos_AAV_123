from scipy.integrate import quad
from numpy import deg2rad, rad2deg, linspace
from math import *
import matplotlib.pyplot as plt

RL = 2.7
L = 10

### Funciones

def pre_functions(x,beta):
    R = (RL * (x/L)**(1/3))
    dR = (RL/L**(1/3) * 1/(3*x**(2/3)))
    cos2phi =  ( (-cos(alpha)*dR+sin(alpha)*sin(beta)) * (1+dR**2)**(-0.5) )
    nx = (-dR * (1+dR**2)**-0.5)
    ny = (cos(beta) * (1+dR**2)**-0.5)
    nz = (sin(beta) * (1+dR**2)**-0.5)
    return cos2phi, nx, ny, nz, R

def beta0(x):
    try:
        val = RL/L**(1/3) * (3*tan(alpha)*x**(2/3))**-1
    except:
        val =2*pi
    if val >=1 : return 2*pi
    else: return asin(val)

def funx(beta,x):
    cos2phi, nx, ny, nz, R = pre_functions(x,beta)
    val = (cos2phi * nx * R)
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

def funy(beta,x):
    cos2phi, nx, ny, nz, R = pre_functions(x,beta)
    val = (cos2phi * ny * R)
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

def funz(beta,x):
    cos2phi, nx, ny, nz, R = pre_functions(x,beta)
    val = (cos2phi * nz * R)
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

def funM(beta,x):
    cos2phi, nx, ny, nz, R = pre_functions(x,beta)
    val = (cos2phi * nz * R * x)
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

#Integracion

#Cx
def cx1_int(x):
    val = quad(funx, 0, beta0(x), args=(x))[0]
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

def cx2_int(x):
    val = quad(funx, pi-beta0(x), 2*pi, args=(x))[0]
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

# Cy
def cy1_int(x):
    val = quad(funy, 0, beta0(x), args=(x))[0]
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

def cy2_int(x):
    val = quad(funy, pi-beta0(x), 2*pi, args=(x))[0]
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

# Cz
def cz1_int(x):
    val = quad(funz, 0, beta0(x), args=(x))[0]
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

def cz2_int(x):
    val = quad(funz, pi-beta0(x), 2*pi, args=(x))[0]
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

# Cm0
def cm01_int(x):
    val = quad(funM, 0, beta0(x), args=(x))[0]
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

def cm02_int(x):
    val = quad(funM, pi-beta0(x), 2*pi, args=(x))[0]
    if isinstance(val, complex): print('*** Error in ***',val,'x =', x,'beta =',beta)
    return val

###############################

K = 2

alpha = deg2rad(10)
print('alpha = ', rad2deg(alpha), 'deg')

cx1 = quad(cx1_int, 0, L)[0]
cx2 = quad(cx2_int, 0, L)[0]
cx = (cx1+cx2)*K/(pi*RL**2)
print('Cx =', cx)

cy1 = quad(cy1_int, 0, L)[0]
cy2 = quad(cy2_int, 0, L)[0]
cy = (cy1+cy2)*K/(pi*RL**2)
print('Cy =',cy)

cz1 = quad(cz1_int, 0.1, L)[0]
cz2 = quad(cz2_int, 0.1, L)[0]
cz = (cz1+cz2)*K/(pi*RL**2)
print('Cz =',cz)

cm01 = quad(cm01_int, 0.1, L)[0]
cm02 = quad(cm02_int, 0.1, L)[0]
cm0 = (cm01+cm02)*K/(pi*RL**2)
print('Cm0 =',cm0)

### Plot

a = deg2rad(linspace(0,30,50))
cx_list = []
cz_list = []
cm0_list = []
xcp_list = []

for alpha in a:

    cx1 = quad(cx1_int, 0, L)[0]
    cx2 = quad(cx2_int, 0, L)[0]
    cx = (cx1+cx2)*K/(pi*RL**2)
    cx_list.append(cx)

    cz1 = quad(cz1_int, 0.1, L)[0]
    cz2 = quad(cz2_int, 0.1, L)[0]
    cz = (cz1+cz2)*K/(pi*RL**2)
    cz_list.append(cz)

    cm01 = quad(cm01_int, 0.1, L)[0]
    cm02 = quad(cm02_int, 0.1, L)[0]
    cm0 = (cm01+cm02)*K/(pi*RL**2)
    cm0_list.append(cm0)

    xcp_list.append(cm0/cz)


fig = plt.figure()
plt.plot(rad2deg(a),cx_list)
plt.plot(rad2deg(a),cz_list)
plt.plot(rad2deg(a),cm0_list)
plt.plot(rad2deg(a),xcp_list)
plt.xlabel('alpha[deg]', horizontalalignment='right', x=1.0)
plt.legend(['cx','cz','cm0','xcp[m]'])
plt.grid()
plt.savefig('alphas.pdf')
plt.close()
