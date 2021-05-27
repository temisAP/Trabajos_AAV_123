from math import exp

def mars_atm_altitude(h):

# h en metros
# P en Pascales
# rho en kg/m^3
# T en Kelvin

    if h>7000 and h<100e3:
        T = -23.4 - 0.00222*h + 273.15
        P = 0.699 * exp(-9e-5 * h)
    elif h<7000:
        T = -31 - 0.000998*h + 273.15
        P = 0.699 * exp(-9e-5 * h)
    else:
        T = 2.7
        P = 0

    rho = P/(0.1921*T)
    return P, rho, T


def mars_atm_P(z):
	P, rho, T = mars_atm_altitude(z)
	return P

def mars_atm_rho(z):
	P, rho, T = mars_atm_altitude(z)
	return rho

def mars_atm_T(z):
	P, rho, T = mars_atm_altitude(z)
	return T
