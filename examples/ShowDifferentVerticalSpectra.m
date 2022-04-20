% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

lat0=45; lon0=-150; % Ocean Station PAPA
lat0=-39; lon0=10.5; % Agulhas

atlas = VerticalModeAtlas('/Volumes/MoreStorage/Data/VerticalModeAtlas/VerticalModeAtlas-01.nc');

j_star = 3;
H_norm = 1/sum((j_star+(1:3000)).^(-5/2));
H = @(j) H_norm*(j_star+j).^(-5/2);

[Phi_1,Gamma_1,z] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'exact-hydrostatic');

j_star = 20;
H_norm = 1/sum((j_star+(1:3000)).^(-5/2));
H = @(j) H_norm*(j_star+j).^(-5/2);

[Phi_2,Gamma_2] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'exact-hydrostatic');

j_star = 3;
H_norm = 1/sum((j_star+(1:3000)).^(-2));
H = @(j) H_norm*(j_star+j).^(-2);

[Phi_3,Gamma_3,z] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'exact-hydrostatic');

j_star = 3;
H_norm = 1/sum((j_star^2+(1:3000).^2).^(-5/4));
H = @(j) H_norm*(j_star^2+j^2).^(-5/4);

[Phi_4,Gamma_4,z] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'exact-hydrostatic');


% figure
% subplot(1,2,1)
% plot(Phi_1,z,'LineWidth',2), hold on
% plot(Phi_2,z,'LineWidth',2)
% plot(Phi_3,z,'LineWidth',2)
% ylim([min(z) max(z)])
% 
% subplot(1,2,2)
% plot(Gamma_1,z,'LineWidth',2), hold on
% plot(Gamma_2,z,'LineWidth',2)
% plot(Gamma_3,z,'LineWidth',2)
% ylim([min(z) max(z)])

figure
subplot(1,2,1)
plot(Phi_2./Phi_1,z,'LineWidth',2), hold on
plot(Phi_3./Phi_1,z,'LineWidth',2)
plot(Phi_4./Phi_1,z,'LineWidth',2)
plot([1 1],[z(1) z(end)],'k--')
ylim([min(z) max(z)])

subplot(1,2,2)
plot(Gamma_2./Gamma_1,z,'LineWidth',2), hold on
plot(Gamma_3./Gamma_1,z,'LineWidth',2)
plot(Gamma_4./Gamma_1,z,'LineWidth',2)
plot([1 1],[z(1) z(end)],'k--')
ylim([min(z) max(z)])
xlim([0.5 1.5])
legend('(j,m)=(20,2.5)','(j,m)=(3,2)','Matern')