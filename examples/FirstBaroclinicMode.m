% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 2;
LoadFigureDefaults

minLon = 25;

atlas = VerticalModeAtlas('/Users/jearly/Data/VerticalModeAtlas/VerticalModeAtlas-01.nc');
lat = atlas.latitude;
lon = atlas.longitude;
h = atlas.equivalentDepth();

h(isnan(h)) = 0; % 0 really is a good value, since it means no stratification!



region = [min(lon) max(lon) min(lat) max(lat)];
region = [-180 180 min(lat) max(lat)];
addpath('/Users/jearly/Documents/MATLAB/sandwell');
[topo,topo_lat,topo_lon,offset] = readtopo(region,0.25,0.25);
topo(:,1) = topo(:,2);

iLon = find(minLon<lon,1,'first');
lon = circshift(lon,-iLon,1);
lonLabel = lon;
dl = diff(lon);
dl(dl<0)=dl(dl<0)+360;
lon = [lon(1); lon(1) + cumsum(dl)];
h = circshift(h,-iLon,2);
c = sqrt(9.81*h);
[LAT,LON] = ndgrid(lat,lon);

iLon = find(minLon<topo_lon,1,'first');
topo_lon = circshift(topo_lon,-iLon,1);
dl = diff(topo_lon);
dl(dl<0)=dl(dl<0)+360;
topo_lon = [topo_lon(1); topo_lon(1) + cumsum(dl)];
topo = circshift(topo,-iLon,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure!
%
FigureSize = [50 50 figure_width_2col 350*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];


ax1 = axes;
pcolor(topo_lon,topo_lat,topo), shading interp;

ax2 = axes;
levels = 0:0.2:3.4;
[M,contour_handle] = contour(LON,LAT,c(:,:,1),levels,'ShowText','on');
contour_handle.LineWidth = 2;

linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];

cmap = cat(1,[0 0 0],repmat(linspace(0.3,1,1000).',1,3));
colormap(ax1,flip(cmap));
caxis(ax1,[-5 0])
ylim([-70 70])

set( gca, 'FontSize', figure_axis_tick_size);
ax1.YTick = -60:10:60;

xaxis = cat(2,40:20:180,-160:20:20);
xtick = 40:20:max(lon);
% xtick = interp1(lonLabel,lon,xaxis);
labels = cell(length(xaxis),1);
for iLabel=1:length(labels)
labels{iLabel} = sprintf('%d',xaxis(iLabel));
end
ax1.XTick = xtick;
ax1.XTickLabel = labels;

title(ax1,'First baroclinic mode speed (m/s)','FontSize', figure_title_size, 'FontName', figure_font)

print('FirstBaroclinicModeGlobal.png','-dpng','-r300')