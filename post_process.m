%% house keeping
clc
clear
close all
%%load data
data1 = load('pressure_100.00_101_0.31.dat');
%%define variable radius, pressure and surface tension
press1(:,1)=data1(:,1);
r1(:,1)=1./data1(:,2);
surface_tension(:,1)=data1(:,3);
%%plot
plot(r1,'-x');
hold on;

data2 = load('pressure_100.00_101_0.34.dat');
%%define variable radius, pressure and surface tension
press2(:,1)=data2(:,1);
r2(:,1)=1./data2(:,2);
surface_tension(:,1)=data2(:,3);
%%plot
plot(r2,'-v');

data3 = load('pressure_200.00_201_0.31.dat');
%%define variable radius, pressure and surface tension
press3(:,1)=data3(:,1);
r3(:,1)=1./data3(:,2);
surface_tension(:,1)=data3(:,3);
%%plot
plot(r3,'-d');
 
data4 = load('pressure_200.00_201_0.34.dat');
%%define variable radius, pressure and surface tension
press4(:,1)=data4(:,1);
r4(:,1)=1./data4(:,2);
surface_tension(:,1)=data4(:,3);
%%plot
plot(r4,'-v');

data5 = load('pressure_400.00_401_0.31.dat');
%%define variable radius, pressure and surface tension
press5(:,1)=data5(:,1);
r5(:,1)=1./data5(:,2);
surface_tension(:,1)=data5(:,3);
%%plot
plot(r5,'-^');

data6 = load('pressure_400.00_401_0.34.dat');
%%define variable radius, pressure and surface tension
press6(:,1)=data6(:,1);
r6(:,1)=1./data6(:,2);
surface_tension(:,1)=data6(:,3);
%%plot
plot(r6,'-s');


%%R-P
RP1 = load('C:\Users\s338090\OneDrive - Cranfield University\LBM\Runge-kutta-2nd-ode\P-R-infi\R-P_100_0.31.dat');
rad1(:,1)=RP1(:,2);
plot(rad1,'LineWidth',5);
RP2 = load('C:\Users\s338090\OneDrive - Cranfield University\LBM\Runge-kutta-2nd-ode\P-R-infi\R-P_100_0.34.dat');
rad2(:,1)=RP2(:,2);
plot(rad2,'LineWidth',5);
RP3 = load('C:\Users\s338090\OneDrive - Cranfield University\LBM\Runge-kutta-2nd-ode\P-R-infi\R-P_200_0.31.dat');
rad3(:,1)=RP3(:,2);
plot(rad3,'LineWidth',5);
RP4 = load('C:\Users\s338090\OneDrive - Cranfield University\LBM\Runge-kutta-2nd-ode\P-R-infi\R-P_200_0.34.dat');
rad4(:,1)=RP4(:,2);
plot(rad4,'LineWidth',5);
RP5 = load('C:\Users\s338090\OneDrive - Cranfield University\LBM\Runge-kutta-2nd-ode\P-R-infi\R-P_400_0.31.dat');
rad5(:,1)=RP5(:,2);
plot(rad5,'LineWidth',5);
RP6 = load('C:\Users\s338090\OneDrive - Cranfield University\LBM\Runge-kutta-2nd-ode\P-R-infi\R-P_400_0.34.dat');
rad6(:,1)=RP6(:,2);
plot(rad6,'LineWidth',5);
%%set plot
legend('lx=100-mx=101-0.31','lx=100-mx=101-0.34','lx=200-mx=201-0.31','lx=200-mx=201-0.34','lx=400-mx=401-0.31','lx=400-mx=401-0.34', ...
    'R-P-100-0.31','R-P-100-0.34','R-P-200-0.31','R-P-200-0.34','R-P-400-0.31','R-P-400-0.34');
set(gca,'FontWeight','bold','FontSize',18);
title("Evolution of bubble");
xlabel('t');
ylabel('R');