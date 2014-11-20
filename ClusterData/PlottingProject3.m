clear all
% close all
A=importdata('pos.dat');
time=A(:,1);
totEn=A(:,2);
xEarth=A(:,3);
yEarth=A(:,4);
origo=[0 0];
vxEarth=A(:,5);
vyEarth=A(:,6);
axEarth=A(:,7);
ayEarth=A(:,8);

figure(1);
% hold all
plot(xEarth,yEarth,'-b')
hold on
plot(origo,origo,'or','MarkerSize',10,'MarkerFaceColor','r')
plot(xEarth(end),yEarth(end),'xb')
% xlim([0 60])
legend('Earth', 'Sun')
title('Movement of Earth around the Sun')
ylabel('Distance from Sun (AU)')
xlabel('Distance from Sun (AU)')
% xlim([-1.2 1.2])
% ylim([-1.2 1.2])
axis equal
hold off

% figure(2);
% plot(time,vxEarth,'b');
% hold on
% plot(time,vyEarth,'r');
% legend('vxEarth','vyEarth')
% title('Velocity of Earth')
% xlabel('time (yr)')
% ylabel('Speed (AU/yr)')
% hold off
% 
% figure(3);
% plot(time,axEarth,'b');
% hold on
% plot(time,ayEarth,'r');
% legend('axEarth','ayEarth')
% title('Acceleration of Earth')
% xlabel('time (yr)')
% ylabel('Acceleration (AU/yr^2)')
% hold off

figure(4);
plot(time,totEn,'b');
legend('totEn')
title('Total energy of system')
xlabel('time (yr)')
ylabel('Energy')
hold off

% saveas(1, 'Without_repulsion_w0_01', 'png')