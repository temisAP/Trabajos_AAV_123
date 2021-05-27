classdef MISION < handle
    
    properties
        path
    end
    properties (Constant)
        StefanBoltzmann = 5.6686002E-8;
        absoluteTemperature = 273.15;
        
    end
    methods
        % class constructor
        function obj = ESATAN_TMD_PostProcessing(path)
            % Características del vehículo
            self.beta = beta;    % Coeficiente balístico
            self.E    = E;       % Eficiencia aerodinámica, por defecto es una función constante
            
            % Constantes
            
            self.RT = 6.371e3; %m
            self.g0 = -9.81; %m/s^2
            
            % Gravedad y atmósfera
            % Si no se modifican se usan las predeterminadas
            
            self.rho_fun = self.RHO;
            self.g_fun = self.G;
            
            %
            self.alt_flag = 1
            self.verbose  = 10 %Pasos con output
        end
    end
    
    
end

