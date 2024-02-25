%% System Parameters
m1 = 8; % Mass of the first link in kg
m2 = 8; % Mass of the second link in kg
l1 = 0.5; % Length of the first link in meters
l2 = 1; % Length of the second link in meters
g = 10; % Acceleration due to gravity in m/s^2
%% Compensation Dynamics
syms th2 th1 th1dot th2dot
M=[12+8*cos(th2),8+4*cos(th2);8+4*cos(th2),8];
C=[-4*th2dot*(2*th1dot+th2dot)*sin(th2);4*th1dot^2*sin(th2)];
G=[-80*sin(th1)+sin(th1+th2);-80*sin(th1+th2)];
thdotvec=[th1dot;th2dot];
tau_comp=-inv(M)*(C.*thdotvec+G);
tau_comp=tau_comp(2);
%% Controller Parameters

%k1=1;
%k2=1;
k3=75;
k4=50;

t=out.tout;

out=out.simout;

out=out.data;

plot(t,out(:,1),'b',t,out(:,2),'r');