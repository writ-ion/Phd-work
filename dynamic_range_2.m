clc;
clear all;
close all;

%% Free space path loss

Pt=60000;  %watt
P_dbm=10*log10(Pt)+30;
R_dbm=-123.97;
% c=3*10^8;
f=100; %MHz
% lambda=c/f;
d=1:1:100; %km

FSPL=32.45+24*log10(d)+24*log10(f);
R2=P_dbm+0-FSPL; %no gains considered
plot(d,FSPL,'b-x');
grid on;
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
title('Path loss vs distance');
legend('FSPL');