%%%%%%%%%%%%%%%%%%%%%%%%
%
% Algorithm is as follows:
%   1. Extend the lower boundary to the full depth using the Sandwell
%      topography.
%   2. Manually create knot points switching to a higher number of
%      degrees-of-freedom in the ocean abyss (<2000m).
%   3. Let errors be inversely proportional to the amplitude, N2
%   4. Force f''=0 at the boundaries (in log space).

lat0 = 31; lon0 = -59;
lat0 = -29; lon0 = 1; % few observations
lat0 = 11; lon0 = 142; % deep!
lat0 = 11+22.4/60; lon0 = 142+35.5/60; % stupid deep
lat0 = 38; lon0 = -65.25; % lots of observations
% lat0 = lat0+1; lon0 = lon0+1;

% We extract z_rho just to get the true bounds of the measurements.
[rho,z_rho] = MeanDensityProfileFromLatLon(lat0,lon0,DensityMethod.rho);
[N2,z,rho0] = MeanDensityProfileFromLatLon(lat0,lon0,DensityMethod.stableN2);
N2 = flip(N2);
z = flip(z);

sandwellfile = '/Users/jearly/Documents/MATLAB/sandwell/sandwell.nc';
maxdepth = min(OceanDepthFromLatLon(lat0,lon0,sandwellfile),min(z_rho));
zLim = [maxdepth; max(z_rho)];

K = 3;
constraints = [];

% Group the abyssal (<2500m) ocean measurements together
DF = 2;
abyssalDF = 5;
abyssalDepth = -2000;
abyssalIndex = find(z<abyssalDepth,1,'last');
lastIndex = 1;
if ~isempty(abyssalIndex) && (abyssalIndex-lastIndex)>DF
    % Makes sure there are a minimum of abyssalDF points between knot points
    if (abyssalIndex-lastIndex)/abyssalDF >= 2
        stride=floor((abyssalIndex-lastIndex)/floor((abyssalIndex-1)/abyssalDF));
        midIndex = abyssalIndex + floor((lastIndex - abyssalIndex)/2);
        z_knot = [min(zLim); z(1+stride:stride:abyssalIndex); z(abyssalIndex+DF:DF:end-DF); max(zLim)];
    else
        z_knot = [min(zLim); z(abyssalIndex); z(abyssalIndex+DF:DF:end-DF); max(zLim)];
    end
    % Force the top and bottom to only allow for linear splines (no second derivative).
    if K > 2
        constraints = struct('t',sz([z_knot(1)+(z_knot(2)-z_knot(1))/2; z(end-1)+(z(end)-z(end-1))/2]),'D',[2;2]);
    end
else
    z_knot = [min(zLim); z(1+DF:DF:end-DF); max(zLim)];
end

% Now go ahead and do the correct boundary knot repeats/removals
z_knot = InterpolatingSpline.KnotPointsForPoints(z_knot,K,1);

% These functions are use to transform to and from log-log space of (z,N2)
sz = @(z) flip(log10(1-z));
zs = @(s) flip(1-10.^s);
gN2 = @(N2) flip(log10(N2));
N2g = @(g) flip(10.^g);

% The error in log space is mostly proportional to the magnitude of the function.
% The factor of 4 includes the sqrt, but also just works well.
sigma = flip(10.^(abs(log10(N2/max(N2)))/4 ));


spline = SmoothingSpline(sz(z),gN2(N2),NormalDistribution(1),'t_knot',sz(z_knot),'K',K,'constraints',constraints,'sigma',sigma);
spline.minimize(@(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorFromGCV);
N2function = @(zz) N2g(spline(sz(zz)));

zDense = linspace(min(zLim),max(zLim),5000).';

figure, plot(N2,z,'LineWidth', 2)
hold on, plot(N2function(zDense),zDense,'Color','black','LineWidth',2), xlog, ylog
% legend('raw','spline','smoothed GCV (custom knots)', 'smoothed GCV')

nModes = 100;
N2function = @(zz) N2g(spline(sz(zz)));
im = InternalModesSpectral(N2function,zLim,z_rho,lat0,'nEVP',512,'N2',1,'rho0',rho0);
im.upperBoundary = UpperBoundary.freeSurface;
z_g = im.GaussQuadraturePointsForModesAtFrequency(nModes+1,0);
im = InternalModesSpectral(N2function,zLim,z_g,lat0,'nEVP',512,'nModes',nModes,'N2',1,'rho0',rho0);
im.upperBoundary = UpperBoundary.freeSurface;
im.normalization = Normalization.uMax;
[F,G,h,k,wMaxRatio,kConstantRatio,omegaConstantRatio] = im.ModesAtFrequency(0,'wMax','kConstant','omegaConstant');

h(1:3)
figure, plot(F(:,2),z_g)

return

spline = SmoothingSpline(z,N2,NormalDistribution(1),'lambda',0);
spline.minimize(@(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorFromCV);

figure
plot(N2,z), hold on
plot(spline(z),z)