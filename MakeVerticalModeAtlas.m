netcdfFile = '/Users/jearly/Data/VerticalModeAtlas/VerticalModeAtlas-5d.nc';
lat = -89:5:89;
lon = -179:5:179;
nModes = 100; % number of *fully* resolved modes, include barotropic
version = 1.0;

j = 0:1:(nModes-1);
nLat = length(lat);
nLon = length(lon);
nZ = nModes+1; % assuming w=0 at the bottom
ncPrecision = 'NC_DOUBLE';

cmode = netcdf.getConstant('CLOBBER');
cmode = bitor(cmode,netcdf.getConstant('SHARE'));
cmode = bitor(cmode,netcdf.getConstant('NETCDF4'));
ncid = netcdf.create(netcdfFile, cmode);

% Define the dimensions
latDimID = netcdf.defDim(ncid, 'latitude', nLat);
lonDimID = netcdf.defDim(ncid, 'longitude', nLon);
modeDimID = netcdf.defDim(ncid, 'mode', nModes);
zDimID = netcdf.defDim(ncid, 'depth', nZ);

% Define coordinate variables
latVarID = netcdf.defVar(ncid, 'latitude', ncPrecision, latDimID);
lonVarID = netcdf.defVar(ncid, 'longitude', ncPrecision, lonDimID);
modeVarID = netcdf.defVar(ncid, 'mode', 'NC_INT', modeDimID);
netcdf.putAtt(ncid,latVarID, 'units', 'degrees_north');
netcdf.putAtt(ncid,lonVarID, 'units', 'degrees_east');
netcdf.putAtt(ncid,modeVarID, 'units', 'mode number');

zVarID = netcdf.defVar(ncid, 'z', ncPrecision, [zDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,zVarID, 'long_name', 'quadrature points');
netcdf.putAtt(ncid,zVarID, 'units', 'm');

rhoVarID = netcdf.defVar(ncid, 'rho', ncPrecision, [zDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,rhoVarID, 'long_name', 'potential density');

N2VarID = netcdf.defVar(ncid, 'N2', ncPrecision, [zDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,N2VarID, 'long_name', 'squared buoyancy frequency');
netcdf.putAtt(ncid,N2VarID, 'units', 's^-2');


FVarID = netcdf.defVar(ncid, 'F', ncPrecision, [zDimID,modeDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,FVarID, 'long_name', '(u,v,p) modes');
netcdf.putAtt(ncid,FVarID, 'normalization', 'uMax');
GVarID = netcdf.defVar(ncid, 'G', ncPrecision, [zDimID,modeDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,GVarID, 'long_name', '(w,eta) modes');
netcdf.putAtt(ncid,GVarID, 'normalization', 'wMax');
hVarID  = netcdf.defVar(ncid, 'h', ncPrecision, [modeDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,hVarID, 'long_name', 'equivalent depth');
netcdf.putAtt(ncid,hVarID, 'units', 'm');

FwMaxNormVarID  = netcdf.defVar(ncid, 'FwMaxNorm', ncPrecision, [modeDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,FwMaxNormVarID, 'long_name', 'Conversion to wMax norm for F modes');
FkConstantNormVarID  = netcdf.defVar(ncid, 'FkConstantNorm', ncPrecision, [modeDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,FkConstantNormVarID, 'long_name', 'Conversion to kConstant norm for F modes');
FomegaConstantNormVarID  = netcdf.defVar(ncid, 'FomegaConstantNorm', ncPrecision, [modeDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,FomegaConstantNormVarID, 'long_name', 'Conversion to omegaConstant norm for F modes');

GuMaxNormVarID  = netcdf.defVar(ncid, 'GuMaxNorm', ncPrecision, [modeDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,GuMaxNormVarID, 'long_name', 'Conversion to uMax norm for G modes');
GkConstantNormVarID  = netcdf.defVar(ncid, 'GkConstantNorm', ncPrecision, [modeDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,GkConstantNormVarID, 'long_name', 'Conversion to kConstant norm for G modes');
GomegaConstantNormVarID  = netcdf.defVar(ncid, 'GomegaConstantNorm', ncPrecision, [modeDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,GomegaConstantNormVarID, 'long_name', 'Conversion to omegaConstant norm for G modes');

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Attribution', 'Created by Jeffrey J. Early.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'Version', version);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), 'CreationDate', datestr(datetime('now')));

% End definition mode
netcdf.endDef(ncid);

netcdf.putVar(ncid, latVarID, lat);
netcdf.putVar(ncid, lonVarID, lon);
netcdf.putVar(ncid, modeVarID, j)


shouldRescale = 1;

[LAT,LON] = ndgrid(lat,lon);

g = 9.81;
c = nan([size(LAT),3]);

startTime = datetime('now');
totalProfiles = length(lat)*length(lon);
iFinished = 0;
iSkipped = 0;

for iLat = 1:length(lat)
    for iLon = 1:length(lon)
        if mod(iLon,5) == 0
            timePerStep = (datetime('now')-startTime)/(iFinished+0.1*iSkipped); % we assume skipping takes 1/10th the time
            timeRemaining = 0.66*(totalProfiles-(iFinished+iSkipped))*timePerStep; % And then assume 1/3rds of the remaining profiles will be skipped
            fprintf('\tcomputed %d, skipped %d, so %d of %d finished. Estimated finish time %s (%s from now)\n', iFinished, iSkipped, iFinished+iSkipped, totalProfiles, datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
        end

        lat0 = lat(iLat);
        lon0 = lon(iLon);
%         [rho,z_rho] = MeanDensityProfileFromLatLon(lat0,lon0,DensityMethod.rho);
%         [z_rho,indices] = sort(z_rho);
%         rho = rho(indices);
        [N2,z,rho0] = MeanDensityProfileFromLatLon(lat0,lon0,DensityMethod.stableN2);
        if length(N2) > 10
            [z,indices] = sort(z);
            N2 = N2(indices);
            sandwellfile = '/Users/jearly/Documents/MATLAB/sandwell/sandwell.nc';
            % Need to fix this... this isn't a good assumption. It really
            % messes up near boundaries.
            maxdepth = min(OceanDepthFromLatLon(lat0,lon0,sandwellfile),min(z));
            zLim = [maxdepth; 0];
        
            iFinished = iFinished + 1;
            % N2function = @(zz) interp1(z,N2,zz,'linear','extrap');
            N2function = SmoothN2(N2,z,zLim);

            nEVP = 512;
            isTooFewPoints = 1;
            while (isTooFewPoints == 1)
                isTooFewPoints = 0;
                im = InternalModesSpectral(N2function,zLim,z,lat0,'nEVP',nEVP,'N2',1,'rho0',rho0);
                im.upperBoundary = UpperBoundary.freeSurface;
                try
                z_g = im.GaussQuadraturePointsForModesAtFrequency(nZ,0);
                catch ME
                    switch ME.identifier
                        case 'GLOceanKit:NeedMorePoints'
                            isTooFewPoints = 1;
                            nEVP = nEVP + 64;
                            warning('Too few points. Adding more.')
                        otherwise
                            rethrow(ME)
                    end
                end
            end
            z_g(end)=0;
            im = InternalModesSpectral(N2function,zLim,z_g,lat0,'nEVP',512,'nModes',nModes,'N2',1,'rho0',rho0);
            im.upperBoundary = UpperBoundary.freeSurface;
            im.normalization = Normalization.uMax;
            [F,G,h,k,wMaxRatio,kConstantRatio,omegaConstantRatio] = im.ModesAtFrequency(0,'wMax','kConstant','omegaConstant');
            
            G = G.*wMaxRatio;

            kappaF = InternalModes.ConditionNumberAsFunctionOfModeNumberForModeIndices(F,1:nModes);
            kappaG = InternalModes.ConditionNumberAsFunctionOfModeNumberForModeIndices(G,1:nModes);
            
            if (kappaF(end)>100 || kappaG(end) > 100)
                fprintf('Bad transformation matrices at (lat,lon)=(%f,%f)',lat0,lon0);
            end

            netcdf.putVar(ncid, zVarID, [0 iLat-1 iLon-1], [nZ 1 1], im.z);
            netcdf.putVar(ncid, rhoVarID, [0 iLat-1 iLon-1], [nZ 1 1], im.rho);
            netcdf.putVar(ncid, N2VarID, [0 iLat-1 iLon-1], [nZ 1 1], im.N2);

            netcdf.putVar(ncid, FVarID, [0 0 iLat-1 iLon-1], [nZ nModes 1 1], F);
            netcdf.putVar(ncid, GVarID, [0 0 iLat-1 iLon-1], [nZ nModes 1 1], G);
            netcdf.putVar(ncid, hVarID, [0 iLat-1 iLon-1], [nModes 1 1], h);

            netcdf.putVar(ncid, FwMaxNormVarID, [0 iLat-1 iLon-1], [nModes 1 1], wMaxRatio);
            netcdf.putVar(ncid, FkConstantNormVarID, [0 iLat-1 iLon-1], [nModes 1 1], kConstantRatio);
            netcdf.putVar(ncid, FomegaConstantNormVarID, [0 iLat-1 iLon-1], [nModes 1 1], omegaConstantRatio);

            netcdf.putVar(ncid, GuMaxNormVarID, [0 iLat-1 iLon-1], [nModes 1 1], 1./wMaxRatio);
            netcdf.putVar(ncid, GkConstantNormVarID, [0 iLat-1 iLon-1], [nModes 1 1], kConstantRatio./wMaxRatio);
            netcdf.putVar(ncid, GomegaConstantNormVarID, [0 iLat-1 iLon-1], [nModes 1 1], omegaConstantRatio./wMaxRatio);
        else
            iSkipped = iSkipped+1;
        end
    end
end

netcdf.close(ncid)
