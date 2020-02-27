%==========================================================
% create the initial condition for the Held-Suarez model
%==========================================================

filename   = 'atmos_average.nc';
filename_o = 'init_cond.nc';

lon   = ncread(filename,'lon');
lat   = ncread(filename,'lat');
pfull = ncread(filename,'pfull');

ucomp = ncread(filename, 'ucomp');
vcomp = ncread(filename, 'vcomp');
temp  = ncread(filename, 'temp' );
ps    = ncread(filename, 'ps'   );

%======================
% take the zonal mean
num_lon   = size(ucomp,1);
num_level = size(ucomp,3);
for k = 1:num_level
  um(:,:,k)=ones(num_lon,1)*squeeze(mean(ucomp(:,:,k),1));
  vm(:,:,k)=ones(num_lon,1)*squeeze(mean(vcomp(:,:,k),1));
  tm(:,:,k)=ones(num_lon,1)*squeeze(mean(temp(:,:,k),1));
end
psm(:,:)=ones(num_lon,1)*squeeze(mean(ps(:,:),1));

%======================
% SAVE DATA as NETCDF
    %Open the file
    ncid = netcdf.create(filename_o,'NC_WRITE');
 
        %Define the dimensions
        dimidt     = netcdf.defDim(ncid,'time',netcdf.getConstant('UNLIMITED'));
        dimidpfull = netcdf.defDim(ncid,'pfull',size(um,3));
        dimidlat   = netcdf.defDim(ncid,'lat',  size(um,2));
        dimidlon   = netcdf.defDim(ncid,'lon',  size(um,1));
 
        %Define IDs for the dimension variables (time,pfull,lat,lon)
        time_ID  = netcdf.defVar(ncid, 'time', 'double', [dimidt]);
        pfull_ID = netcdf.defVar(ncid, 'pfull','double', [dimidpfull]);
        lat_ID   = netcdf.defVar(ncid, 'lat',  'double', [dimidlat]);
        lon_ID   = netcdf.defVar(ncid, 'lon',  'double', [dimidlon]);
 
        %Define the main variable ()
        ucomp_ID  = netcdf.defVar(ncid, 'u', 'float', [dimidlon dimidlat dimidpfull dimidt]);
        vcomp_ID  = netcdf.defVar(ncid, 'v', 'float', [dimidlon dimidlat dimidpfull dimidt]);
        temp_ID   = netcdf.defVar(ncid, 't', 'float', [dimidlon dimidlat dimidpfull dimidt]);
        ps_ID     = netcdf.defVar(ncid, 'ps','float', [dimidlon dimidlat dimidt]);

        tracer1_ID  = netcdf.defVar(ncid, 'tracer1', 'float', [dimidlon dimidlat dimidpfull dimidt]);
        tracer2_ID  = netcdf.defVar(ncid, 'tracer2', 'float', [dimidlon dimidlat dimidpfull dimidt]); 

    %We are done defining the NetCdf
    netcdf.endDef(ncid);
 
        %Then store the dimension variables in
        netcdf.putVar(ncid,time_ID, 0, 1, 1);
        netcdf.putVar(ncid,pfull_ID, pfull);
        netcdf.putVar(ncid,lat_ID,   lat);
        netcdf.putVar(ncid,lon_ID,   lon);
 
        %Then store my main variable
	netcdf.putVar(ncid,ucomp_ID, um);
        netcdf.putVar(ncid,vcomp_ID, vm);
        netcdf.putVar(ncid,temp_ID,  tm);
        netcdf.putVar(ncid,ps_ID,    psm);

        netcdf.putVar(ncid,tracer1_ID, zeros(size(um)));
        netcdf.putVar(ncid,tracer2_ID, zeros(size(um)));
 
    %We're done, close the netcdf
    netcdf.close(ncid)

exit
