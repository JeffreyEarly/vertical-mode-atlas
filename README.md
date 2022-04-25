Vertical Mode Atlas
==============

Global ocean atlas of the vertical modes, including equivalent depth and Rossby number

Quick start
------------

The atlas is initialized by pointing to the NetCDF data that houses the vertical modes and stratification. So, on my machine the one-degree atlas is initialized with,
```matlab
atlas = VerticalModeAtlas('/Volumes/Samsung_T5/VerticalModeAtlas/VerticalModeAtlas-01.nc');
```

After that, you can then query the Atlas for data or to make a few plots. For example,
```matlab
    atlas.ShowLowestModes(lat0,lon0)
    atlas.ShowLowestWKBModes(lat0,lon0)
```
will plot the lowest four modes for a given lat/lon from the numerical solutions and the WKB approximated solutions.

The squared buoyancy frequency can also be fetched,
```matlab
    [N2,z] = atlas.N2(lat0,lon0);
```

Various approximations to the vertical structure functions are fetched as follows,
```matlab
j_star = 3;
H_norm = 1/sum((j_star+(1:3000)).^(-5/2));
H = @(j) H_norm*(j_star+j).^(-5/2);

[Phi_hs,Gamma_hs,z] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'exact-hydrostatic');
[Phi_wkb,Gamma_wkb] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'wkb-hydrostatic');
[Phi_igm,Gamma_igm] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'igm');
[Phi_gm,Gamma_gm] = atlas.VerticalStructureFunctionsWithDistribution(lat0,lon0,H,'gm');
```

Note that you can pass any `H(j)` (as a function handle) to set your own vertical spectrum.

