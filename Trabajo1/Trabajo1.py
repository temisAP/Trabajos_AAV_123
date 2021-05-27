from numpy import deg2rad, rad2deg
from scipy import optimize
from math import tan, sin, cos, pi
from isa import *

# Constantes

Ru = 8.314472 # J/(mol K)

print('--- 1 ---')
# 1. Un gas ideal a temperatura T1 = 80K y presión P1 = 1 atm fluye isentrópicamente a un número de Mach
# M1 = 6,9. La relación de calores específcos del gas es γ = 1,4 y su masa molecular es Mm = 29 g/mol.

T1 = 80 #K
P1 = 1*101325 #Pa
M1 = 6.9
gamma = 1.4
Mm = 29e-3 #kg/mol

# Calcular:
print('*** a ***')
# a) La densidad ρ1 y la energía cinética del gas.

R = Ru / Mm

rho1 = P1/(R*T1)
print('rho1 =',round(rho1,3),'kg/m^3')


cv = 1/(gamma-1) * R

e1 = cv * T1
print('e1 =',round(e1,3),'J/kg')

print('*** b ***')
# b) La temperatura y presión de remanso.

a1 = (gamma*R*T1)**0.5
v1 = M1*a1

cp = gamma/(gamma-1) * R

T01 = T1 + v1**2/(2*cp)
print('T01 =',round(T01,3),'K')

P01 = P1*(T01/T1)**(gamma/(gamma-1))
print('P01 =',round(P01,3),'Pa')

print('*** c ***')
# c) La entalpia estática y la de remanso.

h1 = cp*T1
print('h1 =',round(h1,3),'J/kg')

h01 = cp*T01
print('h01 =',round(h01,3),'J/kg')

print('--- 2 ---')
#2. Un gas ideal que fluje a un número de Mach M1 = 8 se encuentra con una onda de choque normal. La
#temperatura es T1 = 80K y la presión P1 = 1 atm. La relación de calores específcos del gas es γ = 1,4 y su
#masa molecular es Mm = 29 g/mol.

M1 = 8
T1 = 80 #K
P1 = 1*101325 #Pa
gamma = 1.4
Mm = 29e-3 #kg/mol

#Calcular:
print('*** a ***')
#a) La relación entre la temperatura estática corriente abajo de la onda, T2, y la temperatura de remanso
#corriente arriba de la onda de choque normal T01.

p2_p1 = 1 + 2*gamma/(gamma+1)*(M1**2-1)
rho2_rho1 = (gamma+1)*M1**2 /(2+(gamma-1)*M1**2)

T2_T1 = p2_p1/rho2_rho1

T2 = T1*T2_T1

R = Ru / Mm

a1 = (gamma*R*T1)**0.5
v1 = M1*a1

cp = gamma/(gamma-1) * R

T01 = T1 + v1**2/(2*cp)

print('T2/T01 =',round(T2/T01,3))


print('*** b ***')
#b) La temperatura de remanso corriente abajo de la onda de choque normal,T02.

M22 = (1+((gamma-1)/2) * M1**2) / (gamma*M1**2-(gamma-1)/2)
M2 = M22**0.5

a2 = (gamma*R*T2)**0.5
v2 = a2 * M2


T02 = T2 + v2**2/(2*cp)
print('T02 =',round(T02,3),'K')

print('*** c ***')
#c) La relación entre la presión estática corriente abajo de la onda de choque, P2, y el doble de la presión
#dinámica corriente arriba de la onda de choque, ρ1U1^2.

P2 = P1 * p2_p1
rho1 = P1/(R*T1)
r = P2/(2*rho1*v1**2)

print(' P2/(2ρ1U1^2) =',round(r,3))


print('--- 3 ---')
#3. Un gas ideal que fluje a un número de Mach M1 = 8 se encuentra con una onda de choque oblicua que deflecta
#las líneas de corriente un ángulo δ = 20º.
# La temperatura es T1 = 80K y la presión P1 = 1 atm. La relación
#de calores específicos del gas es γ = 1,4 y su masa molecular es Mm = 29 g/mol.

M1 = 8
delta = deg2rad(20) #rad
T1 = 80 #K
P1 = 1*101325 #Pa
gamma = 1.4
Mm = 29e-3 # kg/mol

#Calcular:
print('*** a ***')
#a) La componente normal del número de Mach corriente arriba de la onda de choque y la temperatura
#estática del gas post-shock.

def F(beta):
    izq = tan(delta)
    der = 2/tan(beta) * (M1**2*sin(beta)**2-1)/(M1**2*(gamma+cos(2*beta))+2)
    return izq-der

betaW = optimize.newton(F,0.1)
betaS = optimize.newton(F,3.14/2)

M1nW = M1*sin(betaW)

M1nS = M1*sin(betaS)

M1 = M1nW

p2_p1 = 1 + 2*gamma/(gamma+1)*(M1**2-1)
rho2_rho1 = (gamma+1)*M1**2 /(2+(gamma-1)*M1**2)
T2_T1 = p2_p1/rho2_rho1
T2W = T1 * T2_T1


M1 = M1nS

p2_p1 = 1 + 2*gamma/(gamma+1)*(M1**2-1)
rho2_rho1 = (gamma+1)*M1**2 /(2+(gamma-1)*M1**2)
T2_T1 = p2_p1/rho2_rho1
T2S = T1 * T2_T1

print('BetaW =',round(rad2deg(betaW),3),'deg')
print(' Weak: M1n =',round(M1nW,3))
print(' Weak: T2 =',round(T2W,3),'K')
print('BetaS =',round(rad2deg(betaS),3),'deg')
print(' Strong: M1n =',round(M1nS,3))
print(' Strong: T2 =',round(T2S,3),'K')

print('--- 4 ---')
#4. La figura muestra la velocidad y la altitud del módulo Soyuz TMA durante la maniobra de rentrada. Asumiendo el modelo de Atmósfera Estándar Internacional (véase ESDU 77021 y ESDU77022) y que la longitud
#característica del módulo de descenso de la Soyuz es L = 2,2 m, calcular para z = {90; 60; 20} km de altitud:

L = 2.2 #m
altitudes = [90, 60, 20] #km
velocidades = {90:8e3, 60:6e3, 20:0.5e3}
T0 = 288.15
rho0 = 1.225
p0 = 101325
temperaturas = {90:186.9,60:245.5,20:216.7}
presiones = {90:1.812e-6*p0, 60:2.005e-4*p0, 20:5.403e-2*p0}
densidades = {90:2.789e-6*rho0,60:2.354e-4*rho0,20:7.187e-2*rho0}
R = 287
gamma = 1.4

def visco_din(T):
    mu0 = 18.27e-6 #Pa
    S1 = 120 #K
    T0 = 291.15 #K

    val =  mu0 * (T/T0)**(3/2) * (T0+S1)/(T+S1)

    return val

Re_list=[]
M_list=[]
Kn_list=[]

for z in altitudes:
    p = presiones[z]
    T = temperaturas[z]
    rho = densidades[z]
    mu = visco_din(T)
    v = velocidades[z]
    Re = rho*v*L / mu


    a = (gamma*R*T)**0.5
    M = v/a

    lam = mu / p * (pi*R*T/(2*M))**0.5
    Kn = lam/L

    Re_list.append(Re)
    M_list.append(M)
    Kn_list.append(Kn)

print('*** a ***')
#a) El número de Reynolds de la corriente libre.

for i in range(len(Re_list)): print(altitudes[i],'km : Re =','%.4g' % Re_list[i])

print('*** b ***')
#b) El número de Mach de la corriente libre.

for i in range(len(M_list)): print(altitudes[i],'km : M =','%.4g' % M_list[i])

print('*** c ***')
#c) El número de Knudsen de la corriente libre.

for i in range(len(Kn_list)): print(altitudes[i],'km : Kn =','%.4g' % Kn_list[i])
