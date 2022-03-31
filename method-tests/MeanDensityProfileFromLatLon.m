function [rhoOrN2, z, rho0] = MeanDensityProfileFromLatLon(lat0,lon0,method)
% WOA spacing
% 0-100m: 5m spacing
% 100-500m: 25m spacing
% 500-2000m: 50m spacing
% 2000-5500m: 100m spacing
temp_file = 'support/world-ocean-atlas/woa18_decav_t00_04.nc';
salinity_file = 'support/world-ocean-atlas/woa18_decav_s00_04.nc';

lat = double(ncread(temp_file,'lat_bnds'));
lon = double(ncread(temp_file,'lon_bnds'));
depth = double(ncread(temp_file,'depth'));
z = -depth;

latIndex = find(lat0 >= lat(1,:) & lat0 < lat(2,:));
lonIndex = find(lon0 >= lon(1,:) & lon0 < lon(2,:));

temperature = double(squeeze(ncread(temp_file,'t_an',[lonIndex latIndex 1 1],[1 1 Inf Inf])));
SP = double(squeeze(ncread(salinity_file,'s_an',[lonIndex latIndex 1 1],[1 1 Inf Inf])));
% temperature_dd = double(squeeze(ncread(temp_file,'t_dd',[lonIndex latIndex 1 1],[1 1 Inf Inf])));

z(isnan(temperature)) = [];
SP(isnan(temperature)) = [];
temperature(isnan(temperature)) = [];

if isempty(z) || length(z) < 3
    rhoOrN2 = [];
    rho0 = [];
    return;
end

addpath(genpath('support/gsw_matlab_v3_06_11'))

% compute pressure from z
p = gsw_p_from_z(z,lat0);

% compute Absolute Salinity from Practical Salinity
SA = gsw_SA_from_SP(SP,p,lon0,lat0);

switch(method)
    case {DensityMethod.rhoStable,DensityMethod.rhoStableFromN2,DensityMethod.stableN2}
        SA = gsw_stabilise_SA_const_t(SA,temperature,p,lon0,lat0);
    otherwise
end

% Convert in-situ temperature, t, into Conservative Temperature, CT
CT = gsw_CT_from_t(SA,temperature,p);

% Now compute potential density
p_ref = 0;
rho = gsw_rho(SA,CT,p_ref);
rho0 = double(rho(1));

if length(SA)<=10
    method = DensityMethod.rho;
end

switch(method)
    case {DensityMethod.rho, DensityMethod.rhoStable}
        % Nothing to do
        rhoOrN2 = rho;
    case {DensityMethod.rhoFromN2,DensityMethod.rhoStableFromN2}
        [N2,p_mid] = gsw_Nsquared(SA,CT,p);
        z_mid = gsw_z_from_p(p_mid,lat0);
        N2_function = InterpolatingSpline(z_mid,N2);
        g = 9.81; 
        rho_s=(-rho0/g)*cumsum(N2_function);
        rho_s=rho_s + (-rho_s(0)+rho0);
        rhoOrN2 = rho_s(z);
    case {DensityMethod.N2,DensityMethod.stableN2}
        [rhoOrN2,p_mid] = gsw_Nsquared(SA,CT,p);
        z = gsw_z_from_p(p_mid,lat0);
    case DensityMethod.N2function
        [N2,p_mid] = gsw_Nsquared(SA,CT,p);
        z = gsw_z_from_p(p_mid,lat0);
        z_mid = gsw_z_from_p(p_mid,lat0);
        spline = InterpolatingSpline(z_mid,N2);
        rhoOrN2 = @(z) spline(z);
end

end