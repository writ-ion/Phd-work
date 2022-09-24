%% Radar equation & 2-ray model
clc;
close all;
clear all;

%% Input parameters

Pt=60000;   %watt
P_dbm=10*log10(Pt)+30;
R_dbm=-123.93;
Pr=(10^(R_dbm/10))*1000;  %min reception level
Gt=1;
Gr=1;
c=3*10^8;
f=100*10^6;
lambda=c/f;
eta=0.03;
Rt=30000; %Transmitter to target=30 kms
Rr=30000; %Receiver to target=10kms
h_t=200; %transmitter height
h_r=1; %receiver height
d=1000:1000:100000; %distance to receiver

%% Loss for Radar Equation

L1=(Pt*Gt*Gr*(lambda^2)*eta)/(((4*pi)^3)*(Rt^2)*(Rr^2)*Pr);
l1_dB=10*log10(L1);

%% Loss for 2-ray model

L_2ray=40*log10(d)-20*log10(h_t*h_r);
R_2ray=P_dbm+0-L_2ray;
plot(d,L_2ray,'m-*');
grid on;
hold on;

%% Cost 231 Hata model

A=69.55;
B=26.16;
C=39.5; %propagation slope reduced from 44.9 to ensure > 2.5
Cm=-10; %area correction factor open rural area
%small city
% a_hms1=(1.1*log10(f)-0.7)*Rx - (1.56*log10(f)-0.8); 
%large city
a_hms2=3.2*(log10(11.75*h_r))^2-4.97; 
%propagation loss
%small city
% L1=A+B*log10(f)-13.82*log10(Tx)-a_hms1+(C-6.55*log10(Tx))*log10(d)+Cm;
%large city
L=A+B*log10(f/10^6)-13.82*log10(h_t)-a_hms2+(C-6.55*log10(h_t))*log10(d/10^3)+Cm;
%received power
R1=P_dbm+0-L;%no gains considered
plot(d,L,'r:x');

%% Free space path loss

FSPL=32.45+20*log10(d/10^3)+20*log10(f/10^6);
R2=P_dbm+0-FSPL; %no gains considered
plot(d,FSPL,'b-x');
xlabel('Distance');
ylabel('Path Loss (dB)');
title('Path loss vs distance');
legend('Okumura-Hata','FSPL');

%% Received Signal Strength

figure (2);
plot(d,R_2ray,'r-*');
hold on;
grid on;
plot(d,R1,'m:x');
plot(d,R2,'b-x');
xlabel('Distance');
ylabel('Received Signal Strength');
title('Received Signal Strength vs Distance');
legend('2-Ray','Okumura-Hata','FSPL');