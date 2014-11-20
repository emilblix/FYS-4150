clear all
% close all
A=importdata('pos.dat');
time=A(:,1);
totEn=A(:,2);
xSun=A(:,3);
ySun=A(:,4);
origo=[0 0];
xEarth=A(:,6);
yEarth=A(:,7);
xMoon=A(:,9);
yMoon=A(:,10);
xMoon2=A(:,12);
yMoon2=A(:,13);

%%
figure(1);
% hold all
% plot(xSun,ySun,'-or','MarkerSize',10,'MarkerFaceColor','r')
plot(xSun,ySun,'-r')
hold on
plot(xEarth,yEarth,'-b')
plot(xMoon,yMoon,'-g')
% plot(xMoon2,yMoon2,'-r')
% plot(origo,origo,'or','MarkerSize',10,'MarkerFaceColor','r')
% plot(xSun(end),ySun(end),'xb')
% xlim([0 60])
legend('Sun', 'Earth','Moon')
title('Movement of Earth around the Sun')
ylabel('Distance from Sun (AU)')
xlabel('Distance from Sun (AU)')
% xlim([-1.2 1.2])
% ylim([-1.2 1.2])
axis equal
hold off

%%
figure(2)
hold off
p = plot([NaN NaN], [NaN NaN],'r-');
hold on;
h = plot([NaN NaN], [NaN NaN],'b-');
g = plot([NaN NaN], [NaN NaN],'g-');
p1 = plot([NaN NaN], [NaN NaN],'ro','MarkerSize',10,'MarkerFaceColor','r');
h1 = plot([NaN NaN], [NaN NaN],'bo','MarkerSize',5,'MarkerFaceColor','b');
g1 = plot([NaN NaN], [NaN NaN],'go','MarkerSize',5,'MarkerFaceColor','g');

% legend('earth','moon','moon2')

x0 = 2*abs(xEarth(1)-xMoon(1));

for i=1:1:size(xEarth)
    set(p,'XData', xEarth(1:i), 'YData', yEarth(1:i));
    set(h,'XData', xMoon(1:i), 'YData', yMoon(1:i));
    set(g,'XData', xMoon2(1:i), 'YData', yMoon2(1:i));
    set(p1,'XData', xEarth(i), 'YData', yEarth(i));
    set(h1,'XData', xMoon(i), 'YData', yMoon(i));
    set(g1,'XData', xMoon2(i), 'YData', yMoon2(i));
    
    axis([xEarth(i)-x0 xEarth(i)+x0 yEarth(i)-x0 yEarth(i)+x0])
    pause(0.005);
end

%%
% SunDist='Max distance traveled by Sun (AU) = '
% SunDist=max(xSun)-min(xSun)
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