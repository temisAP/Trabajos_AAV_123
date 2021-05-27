RL = 2.7;
L = 10;

alpha = deg2rad(0);

%% Functions 
R = @(x) RL * (x/L)^(1/3);
dR = @(x) RL/L^(1/3) * 1/(3*x^(2/3)); 
cos2phi = @(x,beta) (-cos(alpha)*dR+sin(alpha)*sin(beta)) * (1+dR^2)^-0.5;
n = @(x,beta) [-dR; cos(beta) ; sin(beta)] * (1+dR^2)^-0.5;
nx = @(x,beta) -dR * (1+dR^2)^-0.5;
ny = @(x,beta) cos(beta) * (1+dR^2)^-0.5;
nz = @(x,beta) sin(beta) * (1+dR^2)^-0.5;
beta0 = @(x) asin(RL/L^(1/3) * (3*tan(alpha)*x^(2/3))^-1);

R = @(x) 2.7 * (x/10)^(1/3);
dR = @(x) 2.7/10^(1/3) * 1/(3*x^(2/3)); 
cos2phi = @(x,beta) (-cos(deg2rad(0))*dR+sin(deg2rad(0))*sin(beta)) * (1+dR^2)^-0.5;
n = @(x,beta) [-dR; cos(beta) ; sin(beta)] * (1+dR^2)^-0.5;
nx = @(x,beta) -dR * (1+dR^2)^-0.5;
ny = @(x,beta) cos(beta) * (1+dR^2)^-0.5;
nz = @(x,beta) sin(beta) * (1+dR^2)^-0.5;
beta0 = @(x) asin(2.7/10^(1/3) * (3*tan(deg2rad(0))*x^(2/3))^-1);

%% Integration
fun = @(x,beta) cos2phi * n * R;
funx = @(x,beta) cos2phi * nx * R;
funy = @(x,beta) cos2phi * ny * R;
funz = @(x,beta) cos2phi * nz * R;

cd = integral2(@(x,beta) funx(x,beta),0,L,0,beta0);

%cd = int( int(funx,0,beta0) + int(funx,pi-beta0,2*pi) ,0,L);
