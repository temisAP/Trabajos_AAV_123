close all, clear all,
clc


filename    = 'esfera1.stl';
[p,t,tnorm] = import_stl_fast(filename,1);

figure(2)
patch('vertices',p,'faces',t,'FaceColor','none','edgecolor','k')


U_infty = zeros(size(tnorm));
U_infty(:,2) = 1;
U_infty(:,1) = 0;

Stheta = zeros(length(tnorm),1);

for I=1:length(tnorm)

    Stheta(I,1) = dot(U_infty(I,:),tnorm(I,:))./norm(U_infty(I,:));
    
end

% Cálculo del cp. Método de Newton
cps = Stheta.^2;
cps = [cps;0];

% Cálculo del cp en las zonas de sotavento
I_cpO           = find(Stheta>0);
cps(I_cpO,1)    = 0;

% Cálculo del coeficiente de resistencia (Newton)

cpN = 2*cps;

cd = 0;

% Cálculo del coeficiente de resistencia (Newton-Modificado)

gamma = 1.4;
cpNM = ((gamma+1)/(4*gamma))^(gamma/(gamma+1))*(4/(gamma+1))*cps;

% Cálculo del coeficiente de resistencia (Newton-Busemann) (Esfera)

cpNB = 2*cps + 2/3*(1-cps).^(3/2);

%% Plot cp

% Newton
figure(1)
colormap(jet)
set(gca,'CLim',[0 2])
pat = patch('vertices',p,'faces',t,'FaceColor','b','edgecolor','none')
xlabel('x'),ylabel('y'),zlabel('z'),


set(gca,'CLim',[0 2])
cdata = cpN;
set(pat,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled')
colorbar

% Newton modificado
figure(2)
colormap(jet)
set(gca,'CLim',[0 2])
pat = patch('vertices',p,'faces',t,'FaceColor','b','edgecolor','none')
xlabel('x'),ylabel('y'),zlabel('z'),


set(gca,'CLim',[0 2])
cdata = cpNM;
set(pat,'FaceColor','flat','FaceVertexCData',cdata,'CDataMapping','scaled')
colorbar



