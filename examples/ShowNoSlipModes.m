atlas = VerticalModeAtlas('/Users/jearly/Data/VerticalModeAtlas/VerticalModeAtlas-01.nc');

lat0 = -50.0; lon0 = -25.0;
nQuadPoints = 5;

rho = atlas.rho(lat0,lon0);
[N2,z] = atlas.N2(lat0,lon0);

% Unfortunately, I don't yet have a method for initialization N2 from
% gridded data, so you have to wrap in a function handle. I need to fix this.
im = InternalModesWKBSpectral(@(zin) interp1(z,N2,zin),[min(z) max(z)],z,lat0,'N2',1,'rho0',rho(end));
im.upperBoundary = UpperBoundary.freeSurface;
im.lowerBoundary = LowerBoundary.noSlip;
im.normalization = Normalization.uMax;
zg = im.GaussQuadraturePointsForModesAtFrequency(nQuadPoints,0);
[F,G,h,k,wMax] = im.ModesAtFrequency(0,'wMax');

G = G.*wMax;

figure
subplot(1,3,1)
N2lim = [0; 1.1*max(sqrt(im.N2)*3600/(2*pi))];
plot(repmat(N2lim,1,nQuadPoints),cat(2,zg,zg).','Color','blue','LineWidth', 1), hold on
plot(sqrt(im.N2)*3600/(2*pi),im.z,'Color','black','LineWidth',2)
xlim(N2lim)
xlabel('N (cph)')

subplot(1,3,2)
plot(F(:,1:4),im.z,'LineWidth',2), hold on
xlabel('(u,v,p)-modes')
title(sprintf('Modes at (%.1f,%.1f), h=%.2f, %.2f, %.2f',lat0,lon0,h(1),h(2),h(3)))

subplot(1,3,3)
plot(G(:,1:4),im.z,'LineWidth',2)
xlabel('(w,\eta)-modes')