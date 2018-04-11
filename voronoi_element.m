clear

addpath('iso2mesh');

lambda = 0.4; % wavelength: microns
rmf = 0.15; % packing density
m = 1.8 + i*0.000188; % refractive index
N_elems = 1 % Number of volume elements;
mesh_file = 'mesh_'; % write mesh to file
elem_ka = 10; % volume element sizeparameter
dist = -3; %power law index
pmax = 10*lambda/2/pi/2;
pmin = 10*lambda/2/pi/5; % 



pallo100t
%load pallo200t
%load pallo380t
k=2*pi/lambda; 
coord2 = p'* elem_ka / k;
etopol2=t';

las = 1;

while (las < N_elems+1)
    
    
%power law distribution    
y = rand;
grain_r = ((pmax^(dist+1) - pmin^(dist+1))*y + pmin^(dist+1))^(1/(dist+1));


gr(las) = grain_r;

grain_size = grain_r * k;

d=ones(size(etopol2,2),1);

grain_r
[coord,etopol,d, dp,Npar] = voronoiN2(coord2,etopol2,d,grain_r,rmf);

if(size(etopol,2) > 0) 

    param = ones(1,size(etopol,2))*m.^2;
    filename = strcat(mesh_file, int2str(las), '.h5')

    hdf5write(filename,'/coord',coord);
    hdf5write(filename,'/etopol',int32(etopol),'WriteMode','append');
    hdf5write(filename,'/param_r',real(param),'WriteMode','append' );
    hdf5write(filename,'/param_i',imag(param),'WriteMode','append');

    las = las + 1;
    end 
    
end



% fig = figure;
% figure(1), plotmesh(coord',[etopol',ones(size(dp))])
% %figure(1), plotmesh(coord',[etopol',dp])
% camva('manual');
% camlight('headlight')
% camproj('perspective')
% fig.GraphicsSmoothing = 'on';
% set(fig, 'renderer', 'OpenGL')
% axis vis3d
% axis equal
