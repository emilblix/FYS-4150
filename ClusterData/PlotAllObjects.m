clear all
% close all
A=importdata('pos.dat');
time=A(:,1);
totEn=A(:,2);
nPlanets= (size(A,2)-2)/3

%%
figure(1)
clf(1)
hold all
for i=1:1:nPlanets
    j=3*i;
    X=A(:,j);
    Y=A(:,j+1);
    plot(X,Y);
end
axis equal

%%
figure (2)
clf(2)
plot(time,totEn,'b');
legend('totEn')
title('Total energy of system')
xlabel('time (yr)')
ylabel('Energy')
hold off

%%
dEn=totEn;
dEn(1)=0;
for i=2:1:size(dEn)
   dEn(i)=totEn(i)-totEn(i-1); 
end
 figure(3)
 clf(3)
 dEn(1:3)=[];
 t=time;
 t(1:3)=[];
 plot(t,dEn,'b');
 legend('totEn')
title('Change in total energy of system')
xlabel('time (yr)')
ylabel('Energy')
hold off