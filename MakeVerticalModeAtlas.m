resolution = '04'; % quarter degree
% resolution = '01'; % one degree
stride = 4;

temperature_file = sprintf('support/world-ocean-atlas/woa18_decav_t00_%s.nc',resolution);
salinity_file = sprintf('support/world-ocean-atlas/woa18_decav_s00_%s.nc',resolution);
landsea_file = sprintf('support/world-ocean-atlas/landsea_%s.msk.txt',resolution);

modeAtlasFile = sprintf('/Users/jearly/Data/VerticalModeAtlas/VerticalModeAtlas-%s.nc',resolution);
modeAtlasFile = sprintf('/Volumes/MoreStorage/Data/VerticalModeAtlas/VerticalModeAtlas-noSlip-%s.nc',resolution);

addpath(genpath('support/gsw_matlab_v3_06_11'));

lat = ncread(temperature_file,'lat');
lon = ncread(temperature_file,'lon');
nModes = 100; % number of *fully* resolved modes, include barotropic
version = 1.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Find the appropriate ocean depth
%
depth = ncread(temperature_file,'depth');
depth = cat(1,depth,(depth(end)+100:100:10000).');
landsea = readtable(landsea_file);
landsea = sortrows(landsea,[1,2]);
oceandepth = depth(landsea.Bottom_Standard_Level);
oceandepth = double(reshape(oceandepth,length(lon),length(lat)));

lat = lat(1:stride:end);
lon = lon(1:stride:end);
oceandepth = oceandepth(1:stride:end,1:stride:end);

j = 0:1:(nModes-1);
nLat = length(lat);
nLon = length(lon);
nZ = nModes+1; % assuming w=0 at the bottom

cmode = netcdf.getConstant('CLOBBER');
cmode = bitor(cmode,netcdf.getConstant('SHARE'));
cmode = bitor(cmode,netcdf.getConstant('NETCDF4'));
ncid = netcdf.create(modeAtlasFile, cmode);

ncPrecision = 'NC_FLOAT';
ncPrecision = 'NC_DOUBLE';

% Define the dimensions
boundsDimID = netcdf.defDim(ncid,'nbounds', 2);
latDimID = netcdf.defDim(ncid, 'lat', nLat);
lonDimID = netcdf.defDim(ncid, 'lon', nLon);
modeDimID = netcdf.defDim(ncid, 'mode', nModes);
zDimID = netcdf.defDim(ncid, 'depth', nZ);

% Copy 'crs' variable
ncid_in = netcdf.open(temperature_file, 'NC_NOWRITE');
varID_in = netcdf.inqVarID(ncid_in, 'crs');
[varName, varType, ~, varNAtts] = netcdf.inqVar( ncid_in, varID_in);
crsVarID = netcdf.defVar(ncid, varName, varType, []);
for i=0:varNAtts-1
	attname = netcdf.inqAttName(ncid_in,varID_in,i);	
	netcdf.putAtt(ncid,crsVarID,attname,netcdf.getAtt(ncid_in,varID_in,attname));
end


% Copy 'lat' variable
varID_in = netcdf.inqVarID(ncid_in, 'lat');
[varName, varType, ~, varNAtts] = netcdf.inqVar( ncid_in, varID_in);
latVarID = netcdf.defVar(ncid, varName, varType, latDimID);
for i=0:varNAtts-1
	attname = netcdf.inqAttName(ncid_in,varID_in,i);	
	netcdf.putAtt(ncid,latVarID,attname,netcdf.getAtt(ncid_in,varID_in,attname));
end

% Copy 'lat_bounds' variable
varID_in = netcdf.inqVarID(ncid_in, 'lat_bnds');
[varName, varType, ~, varNAtts] = netcdf.inqVar( ncid_in, varID_in);
latBndsVarID = netcdf.defVar(ncid, varName, varType, [boundsDimID,latDimID]);
for i=0:varNAtts-1
	attname = netcdf.inqAttName(ncid_in,varID_in,i);	
	netcdf.putAtt(ncid,latBndsVarID,attname,netcdf.getAtt(ncid_in,varID_in,attname));
end

% Copy 'lat' variable
varID_in = netcdf.inqVarID(ncid_in, 'lon');
[varName, varType, ~, varNAtts] = netcdf.inqVar( ncid_in, varID_in);
lonVarID = netcdf.defVar(ncid, varName, varType, lonDimID);
for i=0:varNAtts-1
	attname = netcdf.inqAttName(ncid_in,varID_in,i);	
	netcdf.putAtt(ncid,lonVarID,attname,netcdf.getAtt(ncid_in,varID_in,attname));
end

% Copy 'lat_bounds' variable
varID_in = netcdf.inqVarID(ncid_in, 'lon_bnds');
[varName, varType, ~, varNAtts] = netcdf.inqVar( ncid_in, varID_in);
lonBndsVarID = netcdf.defVar(ncid, varName, varType, [boundsDimID,lonDimID]);
for i=0:varNAtts-1
	attname = netcdf.inqAttName(ncid_in,varID_in,i);	
	netcdf.putAtt(ncid,lonBndsVarID,attname,netcdf.getAtt(ncid_in,varID_in,attname));
end


% Create 'mode' variable
modeVarID = netcdf.defVar(ncid, 'mode', 'NC_INT', modeDimID);
netcdf.putAtt(ncid,modeVarID, 'units', 'mode number');

% Create 'z' variable
zVarID = netcdf.defVar(ncid, 'z', ncPrecision, [zDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,zVarID, 'long_name', 'quadrature points');
netcdf.putAtt(ncid,zVarID, 'units', 'm');

% Create 'rho' variable---double precision, because we lose too many sig figs otherwise
rhoVarID = netcdf.defVar(ncid, 'rho', 'NC_DOUBLE', [zDimID,latDimID,lonDimID]);
netcdf.putAtt(ncid,rhoVarID, 'long_name', 'potential density');

% Create 'N2' variable---also double precision, althogh this is likely not
% necesssary.
N2VarID = netcdf.defVar(ncid, 'N2', 'NC_DOUBLE', [zDimID,latDimID,lonDimID]);
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

netcdf.putVar(ncid, crsVarID, ncread(temperature_file,'crs'));
netcdf.putVar(ncid, latVarID, lat);
latBnds = ncread(temperature_file,'lat_bnds');
latBnds = latBnds(:,1:stride:end);
netcdf.putVar(ncid, latBndsVarID, latBnds);
netcdf.putVar(ncid, lonVarID, lon);
lonBnds = ncread(temperature_file,'lon_bnds');
lonBnds = lonBnds(:,1:stride:end);
netcdf.putVar(ncid, lonBndsVarID, lonBnds);
netcdf.putVar(ncid, modeVarID, j)

startTime = datetime('now');
totalProfiles = sum(oceandepth(:)>0);
iFinished = 0;
iSkipped = 0;

for iLat = 1:length(lat) % 136:136 
    for iLon = 1:length(lon) % 31:31
        if mod(iLon,15) == 0
            timePerStep = (datetime('now')-startTime)/(iFinished+1); % we assume skipping takes 1/10th the time
            timeRemaining = (totalProfiles-iFinished)*timePerStep; % And then assume 1/3rds of the remaining profiles will be skipped
            fprintf('\tcomputed %d of %d finished. Estimated finish time %s (%s from now)\n', iFinished, totalProfiles, datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
        end
        if oceandepth(iLon,iLat) == 0
            iSkipped = iSkipped+1;
            continue;
        end
        lat0 = double(lat(iLat));
        lon0 = double(lon(iLon));

        z = -ncread(temperature_file,'depth');
        temperature = squeeze(ncread(temperature_file,'t_an',[(iLon-1)*stride+1 (iLat-1)*stride+1 1 1],[1 1 Inf Inf]));
        SP = squeeze(ncread(salinity_file,'s_an',[(iLon-1)*stride+1 (iLat-1)*stride+1 1 1],[1 1 Inf Inf]));
        
        z(isnan(temperature)) = [];
        SP(isnan(temperature)) = [];
        temperature(isnan(temperature)) = [];

        % Optimization toolbox features require double
        z = double(z);
        temperature = double(temperature);
        SP = double(SP);
        
        if isempty(z) || length(z) < 10
            iSkipped = iSkipped+1;
        else
            iFinished = iFinished + 1;

            % compute pressure from z
            p = gsw_p_from_z(z,lat0);

            % compute Absolute Salinity from Practical Salinity
            SA = gsw_SA_from_SP(SP,p,lon0,lat0);
            SA = gsw_stabilise_SA_const_t(SA,temperature,p,lon0,lat0);

            % Convert in-situ temperature, t, into Conservative Temperature, CT
            CT = gsw_CT_from_t(SA,temperature,p);

            % Now compute potential density
            p_ref = 0;
            rho = gsw_rho(SA,CT,p_ref);
            rho0 = double(rho(1));
            
            % Now compute the local N2.
            [N2,p_mid] = gsw_Nsquared(SA,CT,p);
            z = gsw_z_from_p(p_mid,lat0);

            % Make sure that we are sorted
            [z,indices] = sort(z);
            N2 = N2(indices);

            % What's the actual lower-boundary here?
            maxdepth = min(oceandepth(iLon,iLat),min(z));
            zLim = [maxdepth; 0];
            
            % Compute a smooth N2 profile, all the way to the
            % lower-boundary
            N2function = SmoothN2(N2,z,zLim);

            nEVP = 256;
            maxEVP = 1024;
            isTooFewPoints = 1;
            lastKappaF = zeros(1,nModes);
            lastKappaG = zeros(1,nModes);
            while (isTooFewPoints == 1 && nEVP <= maxEVP)
                isTooFewPoints = 0;
                % In order to use the WKB spectral, we need to turn
                % 'splitting' 'on' in chebfun, otherwise it uses too many
                % points to attempt to deal with our spline fit.
                im = InternalModesWKBSpectral(N2=N2function,zIn=zLim,zOut=z,latitude=lat0, rho0=rho(end), nEVP=nEVP);
                im.upperBoundary = UpperBoundary.freeSurface;
                im.lowerBoundary = LowerBoundary.noSlip;
                try
                    z_g = im.GaussQuadraturePointsForModesAtFrequency(nZ,0);
                catch ME
                    switch ME.identifier
                        case 'GLOceanKit:NeedMorePoints'
                            isTooFewPoints = 1;
                            nEVP = nEVP + 128;
                            fprintf('Quadrature points not resolved, increasing nEVP.\n')
                            continue;
                        otherwise
                            rethrow(ME)
                    end
                end

                z_g(end)=0;
                im = InternalModesWKBSpectral(N2=N2function,zIn=zLim,zOut=z_g,latitude=lat0, rho0=rho(end), nEVP=nEVP, nModes=nModes);
                im.upperBoundary = UpperBoundary.freeSurface;
                im.lowerBoundary = LowerBoundary.noSlip;
                im.normalization = Normalization.uMax;
                [F,G,h,k,wMaxRatio,kConstantRatio,omegaConstantRatio] = im.ModesAtFrequency(0,'wMax','kConstant','omegaConstant');

                G = G.*wMaxRatio;

                kappaF = InternalModes.ConditionNumberAsFunctionOfModeNumberForModeIndices(F,1:nModes);
                kappaG = InternalModes.ConditionNumberAsFunctionOfModeNumberForModeIndices(G,1:nModes);
                
                % We want to try increasing the number of EVP *if* kappa is
                % large *and* our previous attempt at reducing kappa
                % helped.
                maxKappa = max([kappaF;kappaG]);
                maxdKappa = max(max(abs(kappaF-lastKappaF)), max(abs(kappaG-lastKappaG)));
                if (maxKappa > 100 )
                    lastKappaG = kappaG;
                    lastKappaF = kappaF;
                    if maxdKappa < 1
                        break;
                    else
                        isTooFewPoints = 1;
                        nEVP = nEVP + 128;
                        fprintf('Kappa still large, but dKappa is changing, so increasing nEVP.\n')
                        continue;
                    end
                end
            end
            
            if nEVP > maxEVP
                fprintf('Maximum number of points reached, without resolution at (lat,lon)=(%f,%f), kappa=%f, dKappa=%f\n',lat0,lon0,maxKappa,maxdKappa);
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
        end
    end
end

netcdf.close(ncid)
