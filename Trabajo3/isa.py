from math import exp

def isa_altitude(h, F10 = 330, Ap = 0):

# h en metros
# P en Pascales
# rho en kg/m^3
# T en Kelvin

	g = 9.81 #m/s^2
	M = 28.9e-3 #kg/mol
	R = 287
	T0 = 288.15 #K
	P0 = 101325 #Pa
	rho0 = 1.225 #kg/m^3

	if h < 11000:
		landa = -6.5e-3
		T = T0 + landa*h
		P = P0* (T/T0)**(-g/(R*landa))
		rho = rho0 * (T/T0)**(-g/(R*landa)-1)

	elif h < 25000:
		T11 = 216.65
		P11 = 22552
		rho11 = 0.3629
		#landa=0

		T=T11
		P = P11*exp(-g*(h-11000)/(R*T))
		rho=rho11*exp(-g*(h-11000)/(R*T))

	elif h<47000:

		landa=3e-3
		T25=216.65
		P25=2481
		rho25=0.0399

		T=T25+landa*(h-25000)
		P=P25*(T/T25)**(-g/(R*landa))
		rho=rho25*(T/T25)**(-g/(R*landa)-1)

	elif h<180e3 :

		#CIRA model
		a0 =  7.001985e-2
		a1 = -4.336216e-3
		a2 = -5.009831e-3
		a3 =  1.621827e-4
		a4 = -2.471283e-6
		a5 =  1.904383e-8
		a6 = -7.189421e-11
		a7 =  1.060067e-13

		h = h*1e-3

		polyh = ((((((a7*h + a6)*h + a5)*h + a4)*h + a3)*h + a2)*h + a1)*h + a0

		rho = 10**(polyh)

		#Esto no esta bien
		T = 900 + 2.5 * (F10 - 70) + 1.5*Ap
		P= rho*R*T

	elif h<=550e3:

		R = 287

		T = 900 + 2.5 * (F10 - 70) + 1.5*Ap
		mu = 27 - 0.012 * (h - 200)
		H = T / mu
		try:
			rho = 6*10e-10*exp( - (h - 175) / H )
		except:
			rho = 0

		#Esto no esta bien
		P= rho*R*T

	else:
		#print('Atmósfera no apreciable')

		rho=0
		P=0
		T=0

	return P, rho, T

def isa_P(z):
	P, rho, T = isa_altitude(z)
	return P

def isa_rho(z):
	P, rho, T = isa_altitude(z)
	return rho

def isa_T(z):
	P, rho, T = isa_altitude(z)
	return T
