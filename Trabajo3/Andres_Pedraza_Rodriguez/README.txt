mision_class.py
	Contiene la clase "mision"
	Atributos:
		 Las funciones "beta(gamma)" y "E(gamma)"
		 Las constantes relacionadas con el planeta: "RT", "g0", "mu" (por defecto los de la Tierra)
		 Las funciones "rho_fun(z)" y "g_fun(z)" que pueden ser lambda function, funciones definidas en el c�digo o cualquier funci�n externa, siempre y cuando devuelvan un escalar (por defecto son un modelo exponencial de densidad y la gravedad newtoniana)
		 El error admisible "epsilon" (por defecto 0.1)
	M�todos:
		 reentrada([u0,gamma0,r0,theta0],rfin = 6.371e6,tfin = '')
		 	Devuelve: u, gamma, r, theta
			Adem�s guarda como atributos: u, gamma, r, altitud, theta, t
		
		reentrada_analitica([u0,gamma0,r0,theta0])
			Crea y guarda las im�genes de compartiva: U/Ue vs z y n vs z