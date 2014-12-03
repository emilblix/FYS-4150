% Create cluster


N = 100;  % Number of bodies
R0 = 20;  % Radius of cluster in lightyears

mass= normrnd(10,1,N,1);

rvals = 2*rand(N,1)-1;
elevation = asin(rvals);
azimuth = 2*pi*rand(N,1);
radii = R0*(rand(N,1).^(1/3));
[x,y,z] = sph2cart(azimuth,elevation,radii);

M=zeros(N,7);
M(:,1)=mass;
M(:,2)=x;
M(:,3)=y;
M(:,4)=z;

dlmwrite('cluster.txt',M,' ')
% export(M,'file','cluster.txt','Delimiter',' ')

%%
% for i=1:1:N
%     
%     m= normrnd(10,1);
%     
%     %x= -R0 + 2*R0
%     rvals = 2*rand()-1;
%     elevation = asin(rvals);
%     azimuth = 2*pi*rand();
%     radii = R0*(rand().^(1/3));
%     [x,y,z] = sph2cart(azimuth,elevation,radii);
%     
%     
%     
%     
% end
