clear all
close all
A=importdata('Solution.dat');
x=A(:,1);
closedform=A(:,2);
ans=A(:,3);

hold all
plot(x,ans)
plot(x,closedform)
legend('Numeric solution','Exact solution')
hold off
saveas(fig_handle, 'Solution','png')