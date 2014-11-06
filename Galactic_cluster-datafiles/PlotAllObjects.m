clear all
% close all
A=importdata('pos.dat');
time=A(:,1);
totEn=A(:,2);
nPlanets= (size(A,2)-2)/2;

%%
figure(1)
clf(1)
hold all
for i=1:1:nPlanets
    j=2*i+1;
    X=A(:,j);
    Y=A(:,j+1);
    plot(X,Y);
end
axis equal

%%
figure (2)
plot(time,totEn,'b');
legend('totEn')
title('Total energy of system')
xlabel('time (yr)')
ylabel('Energy')
hold off