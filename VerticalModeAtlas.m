classdef VerticalModeAtlas

    properties
        netcdfFile
        latitude
        longitude
    end

    methods
        function self = VerticalModeAtlas(netcdfFile)
            self.netcdfFile = netcdfFile;
            self.latitude = ncread(self.netcdfFile,'lat');
            self.longitude = ncread(self.netcdfFile,'lon');
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
            xlabel('cph')

            subplot(1,3,2)
            plot(F(:,1:4),z,'LineWidth', 2)
            xlabel('(u,v,p) modes')
            title(sprintf('Modes at (%.1f, %.1f), h=%d m, %.2f m, %.2f m, %.2f m',self.latitude(iLat),self.longitude(iLon), round(h(1)),h(2),h(3),h(4) ))

            subplot(1,3,3)
            plot(G(:,1:4),z,'LineWidth', 2)
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

        function [Phi,Gamma,z] = VerticalStructureFunctions(self,lat0,lon0)
            iLat = interp1(self.latitude,1:length(self.latitude),lat0,'nearest');
            iLon = interp1(self.longitude,1:length(self.longitude),lon0,'nearest');
            
            z = ncread(self.netcdfFile,'z',[1 iLat iLon], [Inf 1 1]);
            mode = ncread(self.netcdfFile,'mode');
            F = ncread(self.netcdfFile,'F',[1 1 iLat iLon], [Inf Inf 1 1]);
            G = ncread(self.netcdfFile,'G',[1 1 iLat iLon], [Inf Inf 1 1]);
            h = ncread(self.netcdfFile,'h',[1 iLat iLon], [Inf 1 1]);
            FkConstantNorm = ncread(self.netcdfFile,'FkConstantNorm',[1 iLat iLon], [Inf 1 1]);
            GkConstantNorm = ncread(self.netcdfFile,'GkConstantNorm',[1 iLat iLon], [Inf 1 1]);

            F = F.*(FkConstantNorm.');
            G = G.*(GkConstantNorm.');

            j_star = 3;
            H = (j_star+(1:3000)).^(-5/2);
            H_norm = 1/sum(H);

            Phi = 0*z;
            Gamma = 0*z;
            for iMode = 2:length(mode)
                j = double(mode(iMode));
                Phi = Phi + (1/h(iMode))*F(:,iMode).*F(:,iMode)*H_norm*(j_star+j).^(-5/2); % the 1/h converts it to the const_F_norm
                Gamma = Gamma + G(:,iMode).*G(:,iMode)*H_norm*(j_star+j).^(-5/2);
            end
        end

    end
end