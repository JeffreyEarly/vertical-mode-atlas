classdef VerticalModeAtlas

    properties
        netcdfFile
        latitude
        longitude
    end

    properties (Constant)
        g = 9.81;
        L_gm = 1.3e3; % thermocline exponential scale, meters
        invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
        E_gm = 6.3e-5; % non-dimensional energy parameter
        E = (1.3e3)*(1.3e3)*(1.3e3)*(5.2e-3)*(5.2e-3)*(6.3e-5);
     end

    methods
        function self = VerticalModeAtlas(netcdfFile)
            self.netcdfFile = netcdfFile;
            self.latitude = double(ncread(self.netcdfFile,'lat'));
            self.longitude = double(ncread(self.netcdfFile,'lon'));
        end


        function ShowLowestModes(self,lat0,lon0)
            iLat = interp1(self.latitude,1:length(self.latitude),lat0,'nearest');
            iLon = interp1(self.longitude,1:length(self.longitude),lon0,'nearest');
            
            z = ncread(self.netcdfFile,'z',[1 iLat iLon], [Inf 1 1]);
            N2 = ncread(self.netcdfFile,'N2',[1 iLat iLon], [Inf 1 1]);

            F = ncread(self.netcdfFile,'F',[1 1 iLat iLon], [Inf Inf 1 1]);
            G = ncread(self.netcdfFile,'G',[1 1 iLat iLon], [Inf Inf 1 1]);
            h = ncread(self.netcdfFile,'h',[1 iLat iLon], [Inf 1 1]);

            figure
            subplot(1,3,1)
            plot(sqrt(N2)*3600/(2*pi),z,'LineWidth', 2, 'Color', 'black')
            ylim([min(z) max(z)])
            xlabel('cph')

            subplot(1,3,2)
            plot(F(:,1:4),z,'LineWidth', 2)
            ylim([min(z) max(z)])
            xlabel('(u,v,p) modes')
            title(sprintf('Modes at (%.1f, %.1f), h=%d m, %.2f m, %.2f m, %.2f m',self.latitude(iLat),self.longitude(iLon), round(h(1)),h(2),h(3),h(4) ))

            subplot(1,3,3)
            plot(G(:,1:4),z,'LineWidth', 2)
            ylim([min(z) max(z)])
            xlabel('(w,\eta) modes')
        end

        function ShowLowestWKBModes(self,lat0,lon0)
            [N2,z] = self.N2(lat0,lon0);
            
            N = sqrt(N2);

            xi = cumtrapz(z,sqrt(N2));
            j = 0:3;
            G = (1./sqrt(N)).*sin(xi*j*pi/xi(end));
            F = (sqrt(N)).*cos(xi*j*pi/xi(end));
            h = ((1/self.g)*(xi(end)./(j*pi)).^2);
            
            Fsign = sign(F(end,:));
            F = Fsign.*F./max(abs(F));
            G = Fsign.*G./max(abs(G));
            
            F(:,1) = ones(size(N2));
            G(:,1) = (z-z(1))/(z(end)-z(1));
            h(1) = abs(z(1));

            figure
            subplot(1,3,1)
            plot(sqrt(N2)*3600/(2*pi),z,'LineWidth', 2, 'Color', 'black')
            ylim([min(z) max(z)])
            xlabel('cph')

            subplot(1,3,2)
            plot(F(:,1:4),z,'LineWidth', 2)
            ylim([min(z) max(z)])
            xlabel('(u,v,p) modes')
            title(sprintf('Modes at (%.1f, %.1f), h=%d m, %.2f m, %.2f m, %.2f m',lat0,lon0, round(h(1)),h(2),h(3),h(4) ))

            subplot(1,3,3)
            plot(G(:,1:4),z,'LineWidth', 2)
            ylim([min(z) max(z)])
            xlabel('(w,\eta) modes')            
        end
        
        function [N2,z] = N2(self,lat0,lon0)
            iLat = interp1(self.latitude,1:length(self.latitude),lat0,'nearest');
            iLon = interp1(self.longitude,1:length(self.longitude),lon0,'nearest');
            
            z = ncread(self.netcdfFile,'z',[1 iLat iLon], [Inf 1 1]);
            N2 = ncread(self.netcdfFile,'N2',[1 iLat iLon], [Inf 1 1]);
        end

        function [rho,z] = rho(self,lat0,lon0)
            iLat = interp1(self.latitude,1:length(self.latitude),lat0,'nearest');
            iLon = interp1(self.longitude,1:length(self.longitude),lon0,'nearest');

            z = ncread(self.netcdfFile,'z',[1 iLat iLon], [Inf 1 1]);
            rho = ncread(self.netcdfFile,'rho',[1 iLat iLon], [Inf 1 1]);
        end

        function h = equivalentDepth(self)
            h = squeeze(ncread(self.netcdfFile,'h',[2 1 1], [1 Inf Inf]));
            h(h>1e4)=nan;
        end
        
        function [Phi,Gamma,z] = VerticalStructureFunctionsWithDistribution(self,lat0,lon0,H,approximation)
            % Optional inputs:
            % approximation — 'exact-hydrostatic', 'wkb-hydrostatic', 'igm', 'gm'.
            % verticalModeSpectrum — H(j)
            iLat = interp1(self.latitude,1:length(self.latitude),lat0,'nearest');
            iLon = interp1(self.longitude,1:length(self.longitude),lon0,'nearest');
%             approximation = 'exact-hydrostatic';

            if strcmp(approximation,'exact-hydrostatic')
                z = ncread(self.netcdfFile,'z',[1 iLat iLon], [Inf 1 1]);
                mode = ncread(self.netcdfFile,'mode');
                F = ncread(self.netcdfFile,'F',[1 1 iLat iLon], [Inf Inf 1 1]);
                G = ncread(self.netcdfFile,'G',[1 1 iLat iLon], [Inf Inf 1 1]);
                h = ncread(self.netcdfFile,'h',[1 iLat iLon], [Inf 1 1]);
                FkConstantNorm = ncread(self.netcdfFile,'FkConstantNorm',[1 iLat iLon], [Inf 1 1]);
                GkConstantNorm = ncread(self.netcdfFile,'GkConstantNorm',[1 iLat iLon], [Inf 1 1]);

                F = F.*(FkConstantNorm.');
                G = G.*(GkConstantNorm.');
            elseif strcmp(approximation,'wkb-hydrostatic')
                [N2,z] = self.N2(lat0,lon0);
                
                N = sqrt(N2);
                xi = cumtrapz(z,sqrt(N2));
                j = 0:100; % add a barotropic mode, although we will ignore it
                h = ((1/self.g)*(xi(end)./(j*pi)).^2);
                G = sqrt(2*self.g/xi(end)).*(1./sqrt(N)).*sin(xi*j*pi/xi(end));
                F = sqrt(2*h/xi(end)).* sqrt(N) .*cos(xi*j*pi/xi(end));

                mode = j;
            end

            Phi = 0*z;
            Gamma = 0*z;
            for iMode = 2:length(mode)
                j = double(mode(iMode));
                Phi = Phi + (1/h(iMode))*F(:,iMode).*F(:,iMode)*H(j); % the 1/h converts it to the const_F_norm
                Gamma = Gamma + G(:,iMode).*G(:,iMode)*H(j);
            end
            Phi = self.L_gm*Phi;
            Gamma = (self.invT_gm*self.invT_gm*self.L_gm/self.g)*Gamma;
        end

        function [Phi,Gamma,z] = VerticalStructureFunctions(self,lat0,lon0)
            % Optional inputs:
            % approximation — 'exact-hydrostatic', 'wkb-hydrostatic', 'igm', 'gm'.
            % verticalModeSpectrum — H(j)
            j_star = 3;
            H_norm = 1/sum((j_star+(1:3000)).^(-5/2));
            H = @(j) H_norm*(j_star+j).^(-5/2);

            [Phi,Gamma,z] = VerticalStructureFunctionsWithDistribution(self,lat0,lon0,H);
        end

    end

    methods  (Access = protected)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Error checking and validation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function isValid = validateOmega(~, omega)
            dOmegaVector = diff(omega);
            if any(dOmegaVector<0)
                error('omega must be strictly monotonically increasing.')
            end
            if max(abs(diff(unique(dOmegaVector)))) > 1e-7
                error('omega must be an evenly spaced grid');
            end
            isValid = 1;
        end

        function isValid = validateSpectrumType(~, omega, spectrumType)
            if (any(omega<0)) && strcmp(approximation,'one-sided')
                error('omega contains negative frequencies, yet you requested a one-sided spectrum. This makes no sense. Try again.');
            end
            isValid = any(validatestring(spectrumType,{'one-sided','two-sided'}));
        end

        function isValid = validateApproximations(~, x)
            isValid = any(validatestring(x,{'exact','wkb', 'wkb-hydrostatic', 'gm'}));
        end

        function isValid = validateZ(self, z)
            isValid = all( z >= min(self.z_in) ) && all(z <= max(self.z_in));
        end

        function [z,approximation] = validateVarianceArguments(self,z,varargin)
            p = inputParser;
            addRequired(p,'z',@(x) self.validateZ(x));
            addOptional(p,'approximation','exact',@(x) self.validateApproximations(x));
            parse(p,z,varargin{:})
            z = p.Results.z;
            approximation = p.Results.approximation;
        end

        function [z,omega,approximation,spectrumType] = validateSpectrumArguments(self,z,omega,varargin)
            if length(varargin) < 2 && all(omega>=0)
                spectrumTypeDefault = 'one-sided';
            else
                spectrumTypeDefault = 'two-sided';
            end

            p = inputParser;
            addRequired(p,'z',@(x) self.validateZ(x));
            addRequired(p,'omega',@(x) self.validateOmega(x));
            addOptional(p,'approximation','exact',@(x) self.validateApproximations(x));
            addOptional(p,'spectrumType',spectrumTypeDefault,@(x) self.validateSpectrumType(omega,x));
            parse(p,z,omega,varargin{:})

            z = p.Results.z;
            omega = p.Results.omega;
            approximation = p.Results.approximation;
            spectrumType = p.Results.spectrumType;
        end
    end

end