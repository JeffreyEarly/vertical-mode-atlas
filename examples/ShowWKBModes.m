% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

lat0=45; lon0=-150; % Ocean Station PAPA
lat0=-39; lon0=10.5; % Agulhas

atlas = VerticalModeAtlas('/Volumes/MoreStorage/Data/VerticalModeAtlas/VerticalModeAtlas-01.nc');

atlas.ShowLowestModes(lat0,lon0)
atlas.ShowLowestWKBModes(lat0,lon0)

j_star = 3;
H_norm = 1/sum((j_star+(1:3000)).^(-5/2));
H = @(j) H_norm*(j_star+j).^(-5/2);

[Phi_hs,Gamma_hs,z] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'exact-hydrostatic');
[Phi_wkb,Gamma_wkb] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'wkb-hydrostatic');
[Phi_igm,Gamma_igm] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'igm');
[Phi_gm,Gamma_gm] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'gm');

figure
subplot(1,2,1)
plot(Phi_wkb./Phi_hs,z,'LineWidth',2), hold on
plot(Phi_igm./Phi_hs,z,'LineWidth',2)
plot(Phi_gm./Phi_hs,z,'LineWidth',2)
plot([1 1],[z(1) z(end)],'k--')
ylim([min(z) max(z)])

subplot(1,2,2)
plot(Gamma_wkb./Gamma_hs,z,'LineWidth',2), hold on
plot(Gamma_igm./Gamma_hs,z,'LineWidth',2)
plot(Gamma_gm./Gamma_hs,z,'LineWidth',2)
plot([1 1],[z(1) z(end)],'k--')
ylim([min(z) max(z)])
xlim([0 3])
legend('wkb','igm','gm')


% figure
% subplot(1,2,1)
% plot(Phi_hs,z,'LineWidth',2), hold on
% plot(Phi_wkb,z,'LineWidth',2)
% plot(Phi_igm,z,'LineWidth',2)
% plot(Phi_gm,z,'LineWidth',2)
% ylim([min(z) max(z)])
% 
% subplot(1,2,2)
% plot(Gamma_hs,z,'LineWidth',2), hold on
% plot(Gamma_wkb,z,'LineWidth',2)
% plot(Gamma_igm,z,'LineWidth',2)
% plot(Gamma_gm,z,'LineWidth',2)
% ylim([min(z) max(z)])
% legend('exact','wkb','igm','gm')