clc;
clear all;
close all;

%% Simulation parameters
f=700; %fm frequency MHz
Tx_pow=60000; %60kW
P_dbm=10*log10(Tx_pow)+30;
Tx_pow_1=1000; %1kW
P1_dbm=10*log10(Tx_pow_1)+30;
Tx=248; %transmitter height meters (Digita)
Sx=1; %sensor height
Rx1=30;
Rx2=80;
d=0:.1:.3; %5 km distance
x=1000:1000:80000;

%link budget formula 
%received power(dBm)=Transmitted power(dBm)+gains(dB)-losses(dB)

%%FSPL
FSPL=32.45+20*log10(d)+20*log10(f);

%% Cost 231 Hata model

A = 69.55;
B = 26.16;
%C = 39.5; %propagation slope reduced from 44.9 to ensure > 2.5
Cm = -10; %area correction factor open rural area
C = 34.12;
%small city
a_hms1=(1.1*log10(f)-0.7)*Sx - (1.56*log10(f)-0.8); 
% a_hms2=(1.1*log10(f)-0.7)*Rx1 - (1.56*log10(f)-0.8);
% a_hms3=(1.1*log10(f)-0.7)*Rx2 - (1.56*log10(f)-0.8);
% a_hms=0;
%large city
% a_hms2=3.2*(log10(11.75*Sx))^2-4.97; 
%propagation loss
%sensor height of 1m
L1=A+B*log10(f)-13.82*log10(Tx)-a_hms1+(C-6.55*log10(Tx))*log10(d)+Cm;
%receiver height 30m
% L2=A+B*log10(f)-13.82*log10(Tx)-a_hms+(C-6.55*log10(Tx))*log10(d)+Cm;
%receiver height 80m
% L3=A+B*log10(f)-13.82*log10(Tx)-a_hms+(C-6.55*log10(Tx))*log10(d)+Cm;
%received power
R1=P_dbm+0-L1;%no gains considered
plot(d,R1,'-*');
hold on;
grid on;
% R2=P_dbm+0-L2;
% plot(d,R2,'-*');
% R3=P_dbm+0-L3;
% plot(d,R3,'-*');
R4=P_dbm+0-FSPL;
plot(d,R4,'-*');
xlabel('Distance (km)')
ylabel('Received signal strength (dBm)')
title('Signal strength vs distance')
legend('Okumura-Hata','FSPL');

%% Free space path loss

% FSPL=32.45+20*log10(d)+20*log10(f);
% R4=P_dbm+0-FSPL; %no gains considered
% plot(d,R4);
% xlabel('Distance (km)')
% ylabel('Received signal strength (dBm)')
% title('Signal strength vs distance')
% grid on;

%% 2-ray model
% 
% L_2ray=40*log10(x)+20*log10(f/40)-20*log10(Tx*Sx);
% R_2ray=P_dbm+0-L_2ray;

%% 1kW Tx power Okumura-Hata

% R3=P1_dbm+0-L;
% plot(d,R3,'r:*');
% xlabel('Distance (km)')
% ylabel('Received signal strength (dBm)')
% title('Signal strength vs distance')
% grid on;
% 
% %% 1kW Tx power FSPL
% 
% R4=P1_dbm+0-FSPL;
% plot(d,R4,'b-*');
% xlabel('Distance (km)')
% ylabel('Received signal strength (dBm)')
% title('Signal strength vs distance')
% grid on;
% legend('Okumura-Hata (sensor)','Okumura-Hata (receiver 30)','Okumura-Hata (receiver 80)','FSPL');

%% Path Loss

figure(2)
plot(d,L1,'-o'); %sensor at 1m
hold on;
grid on;
% plot(d,L2,'-*'); %receiver at 30m
% plot(d,L3,'-*'); %receiver at 80m
plot(d,FSPL,'-*');
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
title('Path loss vs distance');
legend('Okumura-Hata','FSPL');

%% Received signal strength

% rss=[R1;R2];
% save ('rss.mat','rss');