clear all
% close all
A=importdata('pos.dat');
%%
time=A(:,1);
kinEn=A(:,2);
potEn=A(:,3);
totEn=kinEn+potEn;
nPlanets= (size(A,2)-3)/3

%%
figure(1)
clf(1)
hold all
for i=1:1:nPlanets
    j=3*i+1;
    X=A(:,j);
    Y=A(:,j+1);
    plot(X,Y);
end
axis equal

%%
figure (2)
clf(2)
plot(time,totEn,'b');
legend('Total energy')
title('Total energy of system')
xlabel('Time')
ylabel('Energy')
hold off

figure (3)
clf(3)
plot(time,kinEn,'b');
legend('Kinetic energy')
title('Kinetic energy energy of system')
xlabel('Time')
ylabel('Energy')
hold off

figure (4)
clf(4)
plot(time,potEn,'b');
legend('Potential energy')
title('Potential energy energy of system')
xlabel('Time')
ylabel('Energy')
hold off
%%
dEn=totEn;
dEn(1)=0;
for i=2:1:size(dEn)
   dEn(i)=totEn(i)-totEn(i-1); 
end
 figure(5)
 clf(5)
 dEn(1:3)=[];
 t=time;
 t(1:3)=[];
 plot(t,dEn,'b');
 legend('\DeltaEnergy')
title('Change in total energy of system')
xlabel('Time')
ylabel('Energy')
hold off