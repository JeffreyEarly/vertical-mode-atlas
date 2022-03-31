atlas = VerticalModeAtlas('/Volumes/MoreStorage/Data/VerticalModeAtlas/VerticalModeAtlas-01.nc');
h = atlas.equivalentDepth();

figure, pcolor(atlas.longitude,atlas.latitude,h), shading flat, colorbar('eastoutside')
