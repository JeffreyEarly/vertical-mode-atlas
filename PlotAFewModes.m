g = 9.81;
M2 = 745.2*60;

shouldRescale = 1;

lat0 = 43.2; lon0 = -125.1;
lat0 = 20; lon0 = -140;
lat0 = 25; lon0 = -130;
lat0 = 31; lon0 = -50;
[rho,z,rho0] = MeanDensityProfileFromLatLon(lat0,lon0,DensityMethod.rho);
if length(rho) < 10
    error('There are fewer than 10 points.');
end
z_out = double(linspace(min(z),max(z),512)');
nPoints = 10;
imSpectral = InternalModesSpectral(double(rho),double(z),z_out,lat0,'nEVP',512);
imSpectral.upperBoundary = UpperBoundary.freeSurface;
imSpectral.normalization = Normalization.kConstant;
z_g = imSpectral.GaussQuadraturePointsForModesAtFrequency(11,0);
z_g2 = imSpectral.GaussQuadraturePointsForModesAtFrequency(101,0);

imSpectral2 = InternalModesSpectral(double(rho),double(z),z_g2,lat0,'nEVP',512);
imSpectral2.upperBoundary = UpperBoundary.freeSurface;
imSpectral2.normalization = Normalization.kConstant;
[F_s,G_s,h_s] = imSpectral2.ModesAtFrequency(2*pi/M2);
return

[F_s,G_s,h_s] = imSpectral.ModesAtFrequency(2*pi/M2);

figure
subplot(1,3,1)
plot(imSpectral.N2,imSpectral.z)
subplot(1,3,2)
plot(G_s(:,1:4),imSpectral.z)
subplot(1,3,3)
plot(F_s(:,1:4),imSpectral.z)

[N2,z,rho0] = MeanDensityProfileFromLatLon(lat0,lon0,DensityMethod.N2function);
if length(rho) < 10
    error('There are fewer than 10 points.');
end
z_out = double(linspace(min(z),max(z),512)');
imSpectral = InternalModes(N2,double([min(z) max(z)]),z_out,lat0,'nEVP',256, 'method', 'spectral', 'N2',1,'rho0',rho0);
[F_s,G_s,h_s] = imSpectral.ModesAtFrequency(2*pi/M2);

figure
subplot(1,3,1)
plot(imSpectral.N2,imSpectral.z)
subplot(1,3,2)
plot(G_s(:,1:4),imSpectral.z)
subplot(1,3,3)
plot(F_s(:,1:4),imSpectral.z)

return

rho_interp = interp1(double(z),double(rho),z_out);
imFD = InternalModes(rho_interp,z_out,z_out,lat0, 'method', 'finiteDifference','orderOfAccuracy',2);
[F_fd,G_fd,h_fd] = imFD.ModesAtFrequency(0);

iMode = 20;
figure
plot(G_s(:,iMode),z_out), hold on, plot(G_fd(:,iMode),z_out)

N = min(length(h_fd),length(h_s));
idx = 1:N;
figure, plot(abs((h_fd(idx)-h_s(idx))./h_s(idx) ) )