# Matlab samples

This article compiles some useful [matlab](../../archive/linux/matlab.md) samples to use with the MOHID converted [netcdf](../../archive/linux/netcdf.md) files.

## Netcdf MOHID vector plot

Sample code for the Portuguese coast application. It integrates several layers of the horizontal velocity. And plots and saves it.

```
file='F:\Results200606\20060601_Hydrodynamic_1.nc';
nc_dump(file)

lon = nc_varget(file,'lon')
lat = nc_varget(file,'lat')
depth = nc_varget(file,'depth')

siz = [length(lon)-1 length(lat)-1];
lon2d = zeros(siz);
lat2d = zeros(siz);

for i=1:siz(1)
for j=1:siz(2)
   lon2d(i,j) = lon(i);
   lat2d(i,j) = lat(j);
end
end

SumVelX = 0;
SumVelY = 0;

for k=35:42

   Vel_xk = nc_varget(file,'Vel_X',[k 1 1],[1 116 176]);
   SumVelX = SumVelX + Vel_xk * (depth(k+1) - depth(k)) / (depth(43) - depth(36));

   Vel_yk = nc_varget(file,'Vel_Y',[k 1 1],[1 116 176]);
   SumVelY = SumVelY + Vel_yk * (depth(k+1) - depth(k)) / (depth(43) - depth(36));

end

%plot arrows
m_proj('miller', 'lon', [min(lon) max(lon)], 'lat', [min(lat) max(lat)]);
m_quiver(lon2d(1:3:end), lat2d(1:3:end), Vel_xk(1:3:end), Vel_yk(1:3:end), 0, 'k');
hold on;
m_grid('box','fancy','linestyle','none','fontsize', 7);
m_usercoast('coastline.mat','patch', [.5 .5 .5]);
%plot scale
m_quiver(-8.5, 37.5, .1,  0.0, 0, 'k');
htv5 = text(.015, .69, '10 cm s^{-1}');
set(htv5,'FontSize',8);

%save figure as a png file
saveas(gcf, 'MyArrowsFigure.png', 'png');
```

## Read and write netcdf files

```
%Read level2 residual fields
%%%%%%%%%%%% 3D results %%%%%%%%%%%%%%%%%%%%%
file = 'BiscayResidual.nc';
u = nc_varget(file,'Vel_X');
v = nc_varget(file,'Vel_Y');
lon = nc_varget(file,'lon');
lat = nc_varget(file, 'lat');
depth = nc_varget(file, 'depth');
level = nc_varget(file, 'Lev_Z');

%Create level1 residual fields
out = 'MOHIDresidual3D.nc';
nc_create_empty(out, nc_write_mode);
nc_add_dimension(out, 'lon', length(lon));
nc_add_dimension(out, 'lat', length(lat));
nc_add_dimension(out, 'deptht', length(depth));

s_ncvars = struct ( 'Name', ...
                       {   'lon',...
                           'lat',...
                           'deptht',...
                           'u',...
                           'v',...
                           'elevation'}, ...
                  'Nctype', ...
                       {   'float',...
                           'float',...
                           'float',...
                           'float',...
                           'float',...
                           'float'}, ...
                  'Dimension', ...
                       {   {'lon'},...
                           {'lat'},...
                           {'deptht'},...
                           {'deptht' 'lat' 'lon'},...
                           {'deptht' 'lat' 'lon'},...
                           {'lat' 'lon'}}, ...
                  'Attribute', ...
                       {   struct( 'Name','units',...
                                   'Value','degrees_east' ), ...
                           struct( 'Name','units',...
                                   'Value','degrees_north' ), ...
                           struct( 'Name','units',...
                                   'Value','m' ), ...
                           struct( 'Name','units',...
                                   'Value','m s-1' ), ...
                           struct( 'Name','units',...
                                   'Value','m s-1' ), ...
                           struct( 'Name','units',...
                                   'Value','m' ), ...
                        })

for i=1:length(s_ncvars)
   nc_addvar(out, s_ncvars(i));
end

nc_dump(out);
```

## See also

- [Matlab](../../archive/linux/matlab.md)
