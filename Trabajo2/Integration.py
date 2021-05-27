from sympy import *
from numpy import deg2rad
from math import pi

RL = 2.7
L = 10

alpha = deg2rad(5)

x, beta = symbols('x beta')

## Funciones

R = Function('R') (RL * (x/L)**(1/3))
dR = Function('dR') (RL/L**(1/3) * 1/(3*x**(2/3)))
cos2phi = Function('cos2phi') ( (-cos(alpha)*dR+sin(alpha)*sin(beta)) * (1+dR**2)**(-0.5) )
nx = Function('nx') (-dR * (1+dR**2)**-0.5)
ny = Function('ny') (cos(beta) * (1+dR**2)**-0.5)
nz = Function('nz') (sin(beta) * (1+dR**2)**-0.5)
beta0 = Function('beta0') (asin(RL/L**(1/3) * (3*tan(alpha)*x**(2/3))**-1))

funx = Function('funx') (cos2phi * nx * R)
funy = Function('funy') (cos2phi * ny * R)
funz = Function('funz') (cos2phi * nz * R)


cx1 = integrate(funx, (beta, 0, beta0), (x, 0, L))
print(cx1)
cx2 = integrate(funx, (beta, pi-beta0, 2*pi), (x, 0, L))
print(cx1)
print('cx =',cx1+cx2)
