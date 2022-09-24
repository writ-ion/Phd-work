clc;
clear all;
close all;

%% Simulation parameters
f=3500; %fm frequency MHz
P_dbm = 62;
d = 0:0.025:10;

%link budget formula 
%received power(dBm)=Transmitted power(dBm)+gains(dB)-losses(dB)

%%FSPL
L=32.45+20*log10(d)+20*log10(f);

figure;
plot(d,L,'-*');
grid on;
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
title('Path loss vs distance');

%% Rx sensitivity
Pr_dbm = P_dbm + 32 - L;
plot(d,Pr_dbm,'-*');
grid on;
xlabel('Distance (km)')
ylabel('Received signal strength (dBm)')
title('Signal strength vs distance')

%% 
% d= 0.350;
% f=1800;
% L=32.45+20*log10(d)+20*log10(f);
