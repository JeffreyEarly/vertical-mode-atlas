function depth = OceanDepthFromLatLon(lat0,lon0,sandwellfile)
lat = ncread(sandwellfile,'lat');
lon = ncread(sandwellfile,'lon');
i = find(lat<lat0,1,'last');
j = find(lon<lon0,1,'last');
depth = 1000*ncread(sandwellfile,'topo',[i j],[1 1]);
end