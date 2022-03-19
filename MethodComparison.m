lat0 = 31; lon0 = -59;
[rho,z,rho0] = MeanDensityProfileFromLatLon(lat0,lon0,DensityMethod.rho);
if length(rho) < 10
    error('There are fewer than 10 points.');
end
z_out = double(linspace(min(z),max(z),512)');
nModes = 100;
nPoints = 101;
im = InternalModesSpectral(double(rho),double(z),z_out,lat0,'nEVP',512);
im.upperBoundary = UpperBoundary.freeSurface;
z_g = im.GaussQuadraturePointsForModesAtFrequency(nPoints,0);

im = InternalModesSpectral(double(rho),double(z),z_g,lat0,'nEVP',512,'nModes',nModes);
im.upperBoundary = UpperBoundary.freeSurface;
im.normalization = Normalization.uMax;
[F_a,G_a,h_a,k_a,wMaxRatio_a,kConstantRatio_a,omegaConstantRatio_a] = im.ModesAtFrequency(0,'wMax','kConstant','omegaConstant');

errorFunction = @(x,y) max(abs(x-y),[],1)./max(abs(y),[],1);
errorTolerance = 1e-2;

methods = cell(1,1);
methods{1} = DensityMethod.rhoStable;
methods{2} = DensityMethod.rhoFromN2;
methods{3} = DensityMethod.rhoStableFromN2;

for iMethod=1:length(methods)
    [rho,z,rho0] = MeanDensityProfileFromLatLon(lat0,lon0,methods{iMethod});
    if length(rho) < 10
        error('There are fewer than 10 points.');
    end
    z_out = double(linspace(min(z),max(z),512)');
    im = InternalModesSpectral(double(rho),double(z),z_out,lat0,'nEVP',512);
    im.upperBoundary = UpperBoundary.freeSurface;
    z_g = im.GaussQuadraturePointsForModesAtFrequency(nPoints,0);

    im = InternalModesSpectral(double(rho),double(z),z_g,lat0,'nEVP',512,'nModes',nModes);
    im.upperBoundary = UpperBoundary.freeSurface;
    im.normalization = Normalization.uMax;
    [F,G,h,k,wMaxRatio,kConstantRatio,omegaConstantRatio] = im.ModesAtFrequency(0,'wMax','kConstant','omegaConstant');

    max_error = max([errorFunction(h,h_a); errorFunction(F,F_a); errorFunction(G,G_a)],[],1);
    max_norm_error = max([errorFunction(wMaxRatio,wMaxRatio_a); errorFunction(kConstantRatio,kConstantRatio_a); errorFunction(omegaConstantRatio,omegaConstantRatio_a)],[],1);
    fprintf('%s has %d modes that match DensityMethod.rho.\n', methods{iMethod}, find(max_error < errorTolerance,1,'last'));
    fprintf('%s has %d norms that match DensityMethod.rho.\n', methods{iMethod}, find(max_norm_error < errorTolerance,1,'last'));
end