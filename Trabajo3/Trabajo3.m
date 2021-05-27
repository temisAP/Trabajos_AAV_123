clc 
clear all
close all

%% Global parameters

global RT; RT = 6.371e6;    %m
global g0; g0 = 9.81;      %m/s^2

%% Mission data

global E;       E = @(gamma) 0.01*gamma; % 20,20,20!
global beta;    beta = @(gamma) 64;

[Y,t] = reentrada(10,[7e3 deg2rad(45) 100e3+RT 0]);



%% Funciones

% Gravedad
function g = G(r)

global g0;
global RT;

z = r-RT;

g = g0 * (RT/(z+RT))^2;

end

% Densidad
function rho = RHO(r)

global RT;
z = r-RT;

if z <= 100e3 %m
    rho0 = 1.2; %kg/m^3
    z0 = 6.7e3; %m
elseif z> 100e3 %m
    rho0 = 3.9e-9; %kg/m^3
    z0 = 60e3; %m
end

try
    rho = rho0 * exp(-z/z0);
catch
    rho = 1000; %kg/m^3
    
end
       
end

% Aquí se pueden añadir otros esquemas
function [Y,T] = reentrada(Dt,Y0)
global RT;

disp('*** Simulación en curso ***')

[Y,T] = RK4(Dt,Y0,RT,'');

end

% RK4
function [Y,T] = RK4(Dt,Y0,rfin,tfin)

Y = Y0;
T = 0;

y = Y0;
t = 0;

k = 0;

while y(3)>= rfin
    
    k = k+1;
    if mod(k,10) == 0
        disp(['Iteracion:',num2str(k),'// Altitud:', num2str(y(3)/1e3), 'km'])
    end
    
    if isa(tfin,'double')
        if t >= tfin
            break
        end
    end
    
    K1 = F(y);
    K2 = F(y+Dt/2 * K1);
    K3 = F(y+Dt/2 * K2);
    K4 = F(y+Dt   * K3);
    
    y1 = y + Dt/6 * (K1+2*K2+2*K3+K4);
    t = t+Dt;
    
    y = y1;
    
    % Por si el último paso hace que se pase
    if y(3)<0
        y(3) = 0;
    end
    
    Y = [Y ; y];
    T = [T ; t];
    
    if k>= 1e3
        disp('Número máximo de pasos alcanzado')
        break
    end
end

end
            
% Ecuaciones de la dinámica de reentrada
function eqns = F(Y)

u = Y(1);
gamma = Y(2);
r = Y(3);
theta = Y(4);

g = G(r);
rho = RHO(r);
try
    beta_val = beta(gamma);
    E_val = E(gamma);
catch
    global beta ; beta_val = beta(gamma);
    global E    ; E_val = E(gamma);
end

eqns  = [-g*sin(gamma)-rho/(2*beta_val) * u^2,...
    1/u * (u^2 * (rho*E_val/(2*beta_val) + cos(gamma)/r) - g*cos(gamma)),...
    u*sin(gamma),...
    1/r * u*cos(gamma)];

end


