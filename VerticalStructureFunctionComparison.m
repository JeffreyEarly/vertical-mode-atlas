% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

lat0=45; lon0=-150;

[Phi,Gamma,z] = atlas.VerticalStructureFunctions(lat0,lon0);

N2 = atlas.N2(lat0,lon0);

L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;
N0 = invT_gm;
g = 9.81;

PhiGM = sqrt(N2)/(L_gm*N0);
GammaGM = g./(L_gm*N0*sqrt(N2));

L = abs(min(z));

PhiExact = Phi./PhiGM;
GammaExact = Gamma./GammaGM;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure!
%
FigureSize = [50 50 figure_width_1col+8 225*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

subplot(1,2,1)
plot([1 1], [-L 0],'--', 'LineWidth', 1.0*scaleFactor, 'Color', 0*[1 1 1])
hold on
plot(PhiExact,z, 'LineWidth', 1.0*scaleFactor, 'Color', 0*[1 1 1])
xlim([0 1.1*max(PhiExact)])
ylabel('depth (m)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
% leg = legend('GM', 'exact');
% leg.Position(1) = 0.3;
% leg.Position(2) = 0.75;
xlim([0 2])
set( gca, 'FontSize', figure_axis_tick_size);
set(gca, 'YTick', 1000*(-5:1:0));
title('$b N_0 N^{-1}(z) \Phi(z)$','Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);


subplot(1,2,2)
plot([1 1], [-L 0],'--', 'LineWidth', 1.0*scaleFactor, 'Color', 0*[1 1 1])
hold on
plot(GammaExact,z, 'LineWidth', 1.0*scaleFactor, 'Color', 0*[1 1 1])
set(gca, 'YTick', []);
xlim([0 1.1*max(GammaExact)])
title('$b N_0 N(z)\Gamma(z)$','Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);

