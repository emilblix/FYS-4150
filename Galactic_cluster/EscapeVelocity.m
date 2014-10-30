clear all
% close all
A=importdata('pos.dat');
origo=[0 0];
time=A(:,1);
totEn=A(:,2);
xSun=A(:,3);
ySun=A(:,4);
xMercury=A(:,5);
yMercury=A(:,6);
xVenus=A(:,7);
yVenus=A(:,8);
xEarth=A(:,9);
yEarth=A(:,10);
xMars=A(:,11);
yMars=A(:,12);
xJupiter=A(:,13);
yJupiter=A(:,14);
xSaturn=A(:,15);
ySaturn=A(:,16);
xUranus=A(:,17);
yUranus=A(:,18);
% xNeptune=A(:,19);
% yNeptune=A(:,20);
% xPluto=A(:,21);
% yPluto=A(:,22);

figure(1);
% hold all
plot(xSun,ySun,'-or','MarkerSize',10,'MarkerFaceColor','r')
hold on
plot(xMercury,yMercury,'-m')
plot(xVenus,yVenus,'-c')
plot(xEarth,yEarth,'-b')
plot(xMars,yMars,'-r')
plot(xJupiter,yJupiter,'-g')
plot(xSaturn,ySaturn,'-m')
plot(xUranus,yUranus,'-c')
% plot(xNeptune,yNeptune,'-b')
% plot(xPluto,yPluto,'-k')

legend('Sun','V0=2*pi','V0=7.0','V0=7.5','V0=8.0','V0=8.5','Analytical solution','V0=9.0')
title('Different start velocities for a planet with start distance 1 AU from the Sun')
ylabel('Distance from Sun (AU)')
xlabel('Distance from Sun (AU)')
% xlim([-1.2 1.2])
% ylim([-1.2 1.2])
axis equal
hold off


