clear all
% close all
A=importdata('3.dat');
origo=[0 0];
time=A(:,1);
totEn=A(:,2)+A(:,3);
xSun=A(:,4);
ySun=A(:,5);
xEarth3=A(:,7);
yEarth3=A(:,8);

A=importdata('4.dat');
xEarth4=A(:,7);
yEarth4=A(:,8);

A=importdata('5.dat');
xEarth5=A(:,7);
yEarth5=A(:,8);

A=importdata('6.dat');
xEarth6=A(:,7);
yEarth6=A(:,8);

A=importdata('7.dat');
xEarth7=A(:,7);
yEarth7=A(:,8);

A=importdata('8.dat');
xEarth8=A(:,7);
yEarth8=A(:,8);

A=importdata('8_5.dat');
xEarth8_5=A(:,7);
yEarth8_5=A(:,8);

A=importdata('8_75.dat');
xEarth8_75=A(:,7);
yEarth8_75=A(:,8);

A=importdata('8_9.dat');
xEarth8_9=A(:,7);
yEarth8_9=A(:,8);

A=importdata('9.dat');
xEarth9=A(:,7);
yEarth9=A(:,8);

figure(1);
% hold all
plot(xSun,ySun,'-or','MarkerSize',10,'MarkerFaceColor','r')
hold on
plot(xEarth4,yEarth4,xEarth5,yEarth5,xEarth6,yEarth6,xEarth7,yEarth7,xEarth8,yEarth8,xEarth8_5,yEarth8_5,xEarth8_75,yEarth8_75,xEarth8_9,yEarth8_9,xEarth9,yEarth9)

legend('Sun','V0=4','V0=5','V0=6','V0=7','V0=8','V0=8.5','V0=8.75','V0=8.9','V0=9')
title('Different start velocities (in AU/yr) for a planet with start distance 1 AU from the Sun')
ylabel('Distance from Sun (AU)')
xlabel('Distance from Sun (AU)')
% xlim([-1.2 1.2])
% ylim([-1.2 1.2])
axis equal
hold off


