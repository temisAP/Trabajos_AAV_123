clear all, close all,clc

% g0 = 9.81;
% rho0 = 1.225;
% zs  = 7.5e3;
% Rt = 6378e3;

colors = {'-r','-g','-b','-m'};

ze = 100e3;
ue = 7.5e3;
gamma = -5*pi/180;


Es = [0 0.25 0.5 1];
legendas = {'\itE \rm = 0' '\itE \rm = 0.25' '\itE \rm = 0.5' '\itE \rm = 1.0'};

% betas = [50 280 500];
betas = [500];

Dt = 0.1;
for I =1:length(betas)
    
    
    beta    = betas(I);
    
    for J = 1:length(Es)
        
        E = Es(J);
        
        yn = [gamma; ue; ze];
        mov = [yn; 0];
        
        while yn(3) > 0
            
            k1 = ecuacioneF(yn, E, beta);
            k2 = ecuacioneF(yn + k1*Dt/2, E, beta);
            k3 = ecuacioneF(yn + k2*Dt/2, E, beta);
            k4 = ecuacioneF(yn + k3*Dt, E, beta);
            
            yn1 = yn + Dt/6*(k1 + 2*k2 + 2*k3 + k4);
            
            aux = [yn1; k1(2)];
            
            mov = [mov aux];
            
            yn = yn1;
            
            %     disp(yn(3)/1000)
            
        end
        mov = mov(:,1:end-1);
        
        t       = (0:length(mov)-1)*Dt;
        gamat   = mov(1,:);
        ut      = mov(2,:);
        zt      = mov(3,:);
        at      = mov(4,:);
        
        matrizSolu = [t' zt' ut' at' gamat'];
        
        entryDinamicSol(J).z0       = ze;
        entryDinamicSol(J).u0       = ue;
        entryDinamicSol(J).gama0    = gamma;
        entryDinamicSol(J).beta0    = beta;
        entryDinamicSol(J).E0       = E;
        entryDinamicSol(J).results  = matrizSolu;
        
        
        figure(1), hold on
        plot(ut/ut(1),zt/1000,colors{J})
        
        
        figure(2), hold on
        plot(gamat*180/pi,zt/1000,colors{J})
        
        figure(3), hold on
        plot(-at/9.81,zt/1000,colors{J})
        
        
        figure(4),
        subplot(3,1,1), hold on
        plot(t/60,gamat*180/pi,colors{J})
        ylabel('\it\gamma')
        subplot(3,1,2), hold on
        plot(t/60,ut/ut(1),colors{J})
        ylabel('\itu/u_e'),
        subplot(3,1,3), hold on
        plot(t/60,zt/1000,colors{J})
        % axis([0 1.05 0 200])
        ylabel('\itz \rm[km]'),
        xlabel('\itt \rm[min]')
        
        
        figure(5),hold on;
        plot(t/60,zt/1000,colors{J})
        
    end
    
    
end



figure(1), hold on
axis([0 1. 0 120])
legend(legendas,'FontSize',6,'Location','northwest'),
ylabel('\itz \rm[km]'),
xlabel('\itU/U_e')
set(gca,'FontSize',10,'FontName','Times new roman','box','on')

figure(2), hold on
% axis([0 1.05 0 120])
legend(legendas,'FontSize',10,'Location','northwest'),
ylabel('\itz \rm[km]'),
xlabel('\it\gamma')
set(gca,'FontSize',10,'FontName','Times new roman','box','on')

figure(3), hold on
legend(legendas,'FontSize',6,'Location','northwest'),
axis([0 16 0 120])
ylabel('\itz \rm[km]'),
xlabel('\itn')
set(gca,'FontSize',10,'FontName','Times new roman','box','on')

figure(4),
subplot(3,1,1), hold on
axis([0 25 -90 5])

ylabel('\it\gamma \rm[º]')
set(gca,'FontSize',10,'FontName','Times new roman','box','on')

subplot(3,1,2), hold on
axis([0 25 0 1])
ylabel('\itu/u_e'),
legend(legendas,'FontSize',10,'Location','northeast'),
set(gca,'FontSize',10,'FontName','Times new roman','box','on')

subplot(3,1,3), hold on
axis([0 25 0 100])
ylabel('\itz \rm[km]'),
xlabel('\itt \rm[min]')
set(gca,'FontSize',10,'FontName','Times new roman','box','on')



figure(5),hold on;
axis([0 25 0 120])
ylabel('\itz \rm[km]'),
xlabel('\itt \rm[min]')
legend(legendas,'FontSize',6,'Location','northeast'),
set(gca,'FontSize',10,'FontName','Times new roman','box','on')

% text(10,90,'\beta = 500','FontSize',10,'FontName','Times new roman')
% text(10,80,'\gamma_e = -5º','FontSize',10,'FontName','Times new roman')

if 0
    grafWidth   =  18;
    grafAR      =  7/18;
    set(5,'PaperUnits','centimeters','PaperSize',[grafWidth grafWidth*grafAR], 'PaperPosition', [0 0 grafWidth grafWidth*grafAR]);
    
    
    rootPath = pwd;
    figPath =   [rootPath(1:end-3) 'fig'];
    
    
    
    
    %     file2save = [ruta filesep  'ejemploFugoide'];
    %     print(5,'-deps2c',file2save)
    %
    file2save	= [figPath filesep 'SolucionNumerica'];
    print(4,'-dpng',file2save)
    
    
    %
    %     rootPath    = pwd;
    %     srcPath     =   [rootPath(1:end-3) 'src'];
    %     file2save   = strcat(srcPath,filesep,'liftingSol_Es');
    %
    %     save(file2save,'entryDinamicSol')
    
end




