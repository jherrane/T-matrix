function [coord,etopol,d,dp,Npar] = voronoiN2(coord,etopol,d, ...
                                              grain_size,rmf)
                           


dp = d;

cx = max(coord(1,:)) + min(coord(1,:));
cy = max(coord(2,:)) + min(coord(2,:));
cz = max(coord(3,:)) + min(coord(3,:));
cc=[cx;cy;cz]*ones(1,size(coord,2));
coord = coord - cc/2;

rad0 = max(max(coord));
coord = coord*1.5;
rad = max(max(coord));

%vol=0;

%for i1 = 1:size(etopol,2)
%    vol = vol + tetra_volume(coord(:,etopol(1,i1)),coord(:,etopol(2,i1)),...
%coord(:,etopol(3,i1)),coord(:,etopol(4,i1)));
                             %end


Nseed = ceil((2*rad)^3 / (4/3*pi*(grain_size)^3));


seeds = (rand(3,Nseed)-0.5)*2*rad;

seeds2 = seeds(:,find(vlen(seeds) < rad));

Nseed = size(seeds2,2);


for i1=1:size(etopol,2)

       cp=(coord(:,etopol(1,i1)) + coord(:,etopol(2,i1)) ...
	     + coord(:,etopol(3,i1)) + coord(:,etopol(4,i1)))/4;

    [m,ind]=min(vlen(cp*ones(1,Nseed)-seeds2));

    d(i1) = ind;
    dp(i1) = ind;
end

node = coord';
elem = etopol';


aa=unique(d);

ra=rand(length(aa),1);
ra2=rand(length(aa),1);

ind= find(ra(d(:)) < (1-sum(rmf)));
ind2 = find(vlen(seeds2(:,d)) > rad0);

ind3=unique([ind;ind2']);

elem(ind3,:)=[];
dp(ind3) = [];

Npar = size(unique(dp),1);

[node,elem]=removeisolatednode(node,elem);

coord=node';
etopol=elem(:,1:4)';

