clc;
clear all;
close all;

%% Parameters

f = 100; %fm frequency MHz
Tx_pow = 60000; %60kW
P_dbm = 10*log10(Tx_pow)+30;

h_t = 248; %transmitter height meters (Digita)
h_r = 0:1:30;

d=10; %km distance from the transmitter

A=69.55;
B=26.16;
C=39.5; %propagation slope
Cm=-10; %area correction factor

%large city
a_hms1=3.2*(log10(11.75*h_r)).^2-4.97; 
L1 = A+B*log10(f)-13.82*log10(h_t)-a_hms1+(C-6.55*log10(h_t))*log10(d)+Cm;
%received power
R1 = P_dbm+0-L1;%no gains considered
figure
plot(L1,h_r,'r-x');
ylabel("Height of sensor (m)")
xlabel("Loss (dB)")
title("Height vs Loss")
grid on;
figure
plot(R1,h_r,'-o');
ylabel("Height of sensor (m)")
xlabel("Received signal strength (dBm)")
title("Height vs Received signal strength")
grid on;

%small city
a_hms2=(1.1*log10(f)-0.7)*h_r - (1.56*log10(f)-0.8);
L2=A+B*log10(f)-13.82*log10(h_t)-a_hms2+(C-6.55*log10(h_t))*log10(d)+Cm;
%received power
R2=P_dbm+0-L2;%no gains considered
figure
plot(L2,h_r,'r-x');
ylabel("Height of sensor (m)")
xlabel("Loss (dB)")
title("Height vs Loss")
grid on;
figure
plot(R2,h_r,'-o');
ylabel("Height of sensor (m)")
xlabel("Received signal strength (dBm)")
title("Height vs Received signal strength")
grid on;