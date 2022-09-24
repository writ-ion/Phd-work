clc;
close all;
clear all;

%%
Pt = 60000;     % TX power in Watt
P_dbm = 10*log10(Pt/(1e-3)); % TX power in dBm
R_dbm = -123.93; % RX sensitivity in dBm
Pr = (1/1000)*(10^(R_dbm/10));  % min reception level in dBm
Gt = 1;         % unity gain of TX = isotropic antenna = no gain
Gr = 1;         % unity gain of RX = isotropic antenna = no gain
c = 3*10^8;     % speed of light
f = 100*10^6;   % utilized frequency 
lambda = c/f;   % wavelength
eta = .88*(lambda^2); % object = dipole (size: lambda/2)
%% Received power

R = ((Pt*Gt*Gr*(lambda^2)*eta) / ( ((4*pi)^3)*Pr))^(1/4);
R_km = R/1000 % distance between TX and target, when signal reflects back to TX
