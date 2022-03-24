lat0 = 31; lon0 = -59;
% lat0 = -29; lon0 = 1; % few observations
% lat0 = 11; lon0 = 142; % deep!
% lat0 = 11+22.4/60; lon0 = 142+35.5/60; % stupid deep
% lat0 = 38; lon0 = -65.25; % lots of observations
[rho,z_rho] = MeanDensityProfileFromLatLon(lat0,lon0,DensityMethod.rho);
[N2,z,rho0] = MeanDensityProfileFromLatLon(lat0,lon0,DensityMethod.stableN2);

sandwellfile = '/Users/jearly/Documents/MATLAB/sandwell/sandwell.nc';

N2 = flip(N2);
z = flip(z);

K = 3;
constraints = [];

sz = @(z) flip(log10(1-z));
zs = @(s) flip(1-10.^s);
gN2 = @(N2) flip(log10(N2));
N2g = @(g) flip(10.^g);

%%%%%%%%%%%%%%%%%%%%%%%%
%
% Algorithm is as follows:
%   1. Extend the lower boundary to the full depth using the Sandwell
%      topography.
%   2. Manually create knot points switching to a higher number of
%      degrees-of-freedom in the ocean abyss (<2000m).
%   3. Let errors be inversely proportional to the amplitude, N2
%   4. Force f''=0 at the boundaries (in log space).

if 1 == 0
    DF = 1;
    z_knot = [min(z_rho); z(2:end-1); max(z_rho)];
    zLim = [min(z_rho); max(z_rho)];
else
    % Group the abyssal (<2500m) ocean measurements together
    DF = 2;
    abyssalDF = 5;
    
    maxdepth = min(OceanDepthFromLatLon(lat0,lon0,sandwellfile),min(z_rho));
    zLim = [maxdepth; max(z_rho)];
%     ztop = z(z>-2500);
%     z_knot = [min(zLim); ztop; max(zLim)];

    abyssalIndex = find(z<-2000,1,'last');
    lastIndex = 1;
    if ~isempty(abyssalIndex) && (abyssalIndex-lastIndex)>DF
        if (abyssalIndex-lastIndex)/abyssalDF >= 2
            stride=floor((abyssalIndex-lastIndex)/floor((abyssalIndex-1)/abyssalDF));
            midIndex = abyssalIndex + floor((lastIndex - abyssalIndex)/2);
            z_knot = [min(zLim); z(1+stride:stride:abyssalIndex); z(abyssalIndex+DF:DF:end-DF); max(zLim)];
        else
            z_knot = [min(zLim); z(abyssalIndex); z(abyssalIndex+DF:DF:end-DF); max(zLim)];
        end
        if K > 2
            constraints = struct('t',sz([z_knot(1)+(z_knot(2)-z_knot(1))/2; z(end-1)+(z(end)-z(end-1))/2]),'D',[2;2]);
        end
    else
        z_knot = [min(zLim); z(1+DF:DF:end-DF); max(zLim)];
    end

    z_knot = InterpolatingSpline.KnotPointsForPoints(z_knot,K,1);

    z_knot_classic = InterpolatingSpline.KnotPointsForPoints([min(zLim); z(2:end-1); max(zLim)],K,DF);

%     z_knot_lowDF= InterpolatingSpline.KnotPointsForPoints([min(zLim); z(2:end-1); max(zLim)],K,DF);
%     z_knot_highDF= InterpolatingSpline.KnotPointsForPoints([min(zLim); z(2:end-1); max(zLim)],K,abyssalDF);
%     abyssalIndex = find(z_knot_highDF<-2000,1,'last');
end

% WOA spacing
% 0-100m: 5m spacing
% 100-500m: 25m spacing
% 500-2000m: 50m spacing
% 2000-5500m: 100m spacing
% zd = flip(-cat(2,0:20:100,350:250:2000,2500:500:5500).');
% zd(zd<maxdepth) = [];
% zd(1) = maxdepth;
% z_knot = InterpolatingSpline.KnotPointsForPoints(zd,K,1);

% z_knot = InterpolatingSpline.KnotPointsForPoints(z_knot,K,DF);
% f = ConstrainedSpline(z,N2,K,z_knot,NormalDistribution(1),struct('global',ShapeConstraint.positive));

f2 = ConstrainedSpline(sz(z),gN2(N2),K,sz(z_knot),NormalDistribution(1),constraints);
% f = @(z) flip(10.^f2(flip(log10(1-z))));
f = @(z) N2g(f2(sz(z)));

zDense = linspace(min(zLim),max(zLim),5000).';
figure, plot(N2,z), hold on, plot(f(zDense),zDense), xlog, ylog, ylim([zDense(1) zDense(end-1)])

sigma = flip(sqrt(10.^(abs(log10(N2/max(N2)))/2 )));

spline = SmoothingSpline(sz(z),gN2(N2),NormalDistribution(1),'t_knot',sz(z_knot),'K',K,'constraints',constraints,'sigma',sigma);
spline.minimize(@(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorFromGCV);
hold on, plot(N2g(spline(sz(zDense))),zDense)
spline = SmoothingSpline(sz(z),gN2(N2),NormalDistribution(1),'K',K,'constraints',constraints,'t_knot',sz(z_knot_classic),'sigma',sigma);
spline.minimize(@(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorFromGCV);
plot(N2g(spline(sz(zDense))),zDense)
legend('raw','spline','smoothed GCV (custom knots)', 'smoothed GCV')

nModes = 100;
N2function = @(zz) f(zz);
N2function = @(zz) N2g(spline(sz(zz)));
im = InternalModesSpectral(N2function,zLim,z_rho,lat0,'nEVP',512,'N2',1,'rho0',rho0);
im.upperBoundary = UpperBoundary.freeSurface;
z_g = im.GaussQuadraturePointsForModesAtFrequency(nModes+1,0);
im = InternalModesSpectral(N2function,zLim,z_g,lat0,'nEVP',512,'nModes',nModes,'N2',1,'rho0',rho0);
im.upperBoundary = UpperBoundary.freeSurface;
im.normalization = Normalization.uMax;
[F,G,h,k,wMaxRatio,kConstantRatio,omegaConstantRatio] = im.ModesAtFrequency(0,'wMax','kConstant','omegaConstant');

h(1:3)
% figure, plot(F(:,50),z_g)

return

spline = SmoothingSpline(z,N2,NormalDistribution(1),'lambda',0);
spline.minimize(@(aSmoothingSpline) aSmoothingSpline.expectedMeanSquareErrorFromCV);

figure
plot(N2,z), hold on
plot(spline(z),z)