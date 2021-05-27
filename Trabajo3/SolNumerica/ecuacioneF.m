function yP = ecuacioneF(yn, E, beta)

g0 = 9.81;
Rt = 6378e3;


gamma   = yn(1);
u       = yn(2);
z       = yn(3);

r = z + Rt;

  
rhoZ = rho(z);

g = g0*(Rt/r)^2;


Lom = rhoZ * u^2 * E /(2*beta);
Dom = rhoZ * u^2 /(2*beta);

gamaP   = 1/u*(Lom  - g*cos(gamma) + u^2/r*cos(gamma));
uP      = -g*sin(gamma) - Dom;
zP      = u*sin(gamma);


yP = [gamaP; uP; zP];



function rhoOut = rho(z)


if z < 150e3
    rho0    = 1.225;
    zs      = 7.524e3;
else
    rho0    = 3.875e-9;
    zs      = 59.06e3;
end

 rhoOut = rho0*exp(-z/zs);
