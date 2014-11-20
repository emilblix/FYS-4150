clear all
% close all
A=importdata('pos.dat');
time=A(:,1);
xSun=A(:,3);
ySun=A(:,4);
origo=[0 0];
xEarth=A(:,5);
yEarth=A(:,6);
xJupiter=A(:,7);
yJupiter=A(:,8);

figure(1);
% hold all
plot(xSun,ySun,'-or','MarkerSize',10,'MarkerFaceColor','r')
hold on
plot(xEarth,yEarth,'-b')
plot(xJupiter,yJupiter,'-g')
% plot(origo,origo,'or','MarkerSize',10,'MarkerFaceColor','r')
% plot(xSun(end),ySun(end),'xb')
% xlim([0 60])
legend('Sun', 'Earth','Jupiter')
title('Movement of Earth around the Sun')
ylabel('Distance from Sun (AU)')
xlabel('Distance from Sun (AU)')
% xlim([-1.2 1.2])
% ylim([-1.2 1.2])
axis equal
hold off

SunDist='Max distance traveled by Sun (AU) = '
SunDist=max(xSun)-min(xSun)
% figure(2);
% plot(time,xEarth,'b');
% hold on
% plot(time,yEarth,'r');
% legend('xEarth','yEarth')
% title('Velocity of Earth')
% xlabel('time (yr)')
% ylabel('Speed (AU/yr)')
% hold off
% 
% figure(3);
% plot(time,xJupiter,'b');
% hold on
% plot(time,yJupiter,'r');
% legend('axEarth','ayEarth')
% title('Acceleration of Earth')
% xlabel('time (yr)')
% ylabel('Acceleration (AU/yr^2)')
% hold off

% saveas(1, 'Without_repulsion_w0_01', 'png')