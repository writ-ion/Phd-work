clc;
clear all;
close all;

R1=1:50; %distance from Tx to sensor in km
R2=1:50; %distance from Tx to sensor in km
R=1:50;
%R=1:50; %distance from sensor to Rx in km
%R=R1+R2; %distance from Tx to Rx
% G1=1; %Transmit antenna gain
% G2=1; %Receiver antenna gain
c=3e8; %in Hertz
f=100e6; %in Hertz
lambda=c/f; %wavelength
% P_tx=60000; %in watts
% P_n=4.003880708000008e-13; %noise in linear scale
sigma=.88*lambda^2; %0.1:0.1:1;
dr = 70; %dynamic range in dB scale
DR = 10^(dr/10);

temp = sqrt(sigma/(4*pi)*DR);

R2 = (sqrt(sigma/(4*pi)*DR).*R')./ R1;
% R = R1'.*R2./(sqrt(sigma/(4*pi)*DR));
mesh(R2)    
xlabel('R1 distance (km)')
ylabel('R2 distance (km)')
zlabel('R distance (km)')

%% new try
syms deltaR
syms R;
c=3e8; %in m/s
f=100e6; %in Hertz
lambda=c/f; %wavelength in m
sigma=.88*lambda^2; %0.1:0.1:1;
dr = 70; %dynamic range in dB scale
DR = 10^(dr/10);

eqn = (deltaR^2)/R + deltaR - sqrt(sigma/4*pi*DR) == 0
fsurf((deltaR^2)/R + deltaR - sqrt(sigma/4*pi*DR))
axis([0 100 0 100 0 100])
solve(eqn,deltaR)