clc;
clear all;
close all;

%% Parameters

h_bs = 2; %height of base station
h_ue = 1.5; %height of user equipment

d = 0.1:0.1:1; % linear distance (in meters)
d_3d = sqrt(d.^2 + (h_bs-h_ue)^2);

f = 60; % IEEE 802.11 ad model

Pt_dbm = 33; % Transmit power in dBm
Pt_lin = (1/1000)*10^(Pt_dbm/10);

shadowing = 1.5;
small_scale_fading = 6; % in dB

%% Calculation of the antenna directivity gain & directivity mismatch

degree = 45;
dir = degree*pi/180;
syms theta;

%Calculation of the directive gain

P_rad = double(int(sin(2*theta),theta,[0 dir])); % Um*cos(theta),int(sin(theta)*cos(theta))
D = 4/P_rad; %antenna directivity
E = 0.6; %antenna efficiency + other losses including polarization mismatch and modulation losses
G = E*D;
G_db = 10*log10(G)

% Calculation of the loss in gain due to directivity (directivity mismatch)

P_rad_loss = double(int((sin(theta)*cos(theta)^2),theta,[dir pi])); % (Um/2)*cos(theta)^2, int(cos(theta)^2*(sin(theta)))
D_loss = 4/P_rad_loss;
L_dB = 10*log10(D_loss)

EIRP = Pt_dbm + G_db % (https://www.60ghz-wireless.com/60ghz-technology/60ghz-band-regulation/) 40 dBm max EIRP

%% Receiver sensitivity

B=100;        %in MHz
B=B*10^6;   %bandwidth in hertz
T=290;      %temperature in kelvin
k=1.38064852*10^-23;     %boltzmann constant

NF=10;      %noise figure in dB    
% SNR=10;     %in dB 
noise_floor = 10*log10(k*T*B) + 30 + NF % (+ SNR) in dBm
noise_linear = (10^(noise_floor/10)); %noise floor in linear scale

%% Propagation model

PL_los_IEEE = 32.5 + 20*log10(f) + 20*log10(d) + small_scale_fading; %no shadow fading mentioned
figure
plot(d,PL_los_IEEE,'-.')
hold on
grid on

PL_nlos_IEEE = 44.2 + 20*log10(f) + 18*log10(d) + shadowing + small_scale_fading; %shadow fading = 1.5
plot(d,PL_nlos_IEEE,'--')
legend('LOS','NLOS')
xlabel('distance (m)')
ylabel('Path loss + Shadowing + Small scale fading')
title ('Loss as a function of distance, IEEE 802.11 ad')
hold off

%% Received Power at the Backscatter device

Pr_nlos_BD = EIRP - PL_nlos_IEEE; % Transmitted EIRP - Path loss in the NLOS environment
Pr_los_BD = EIRP - PL_los_IEEE;
figure
plot(d,Pr_nlos_BD,'-d')
hold on
plot(d,Pr_los_BD,'-s')
xlabel('distance (m)')
ylabel('Received/Reflected signal power (dBm)')
title('Received signal power at the Backscatter device')
legend('NLOS-received power','LOS-received power')
grid on

%% Received power at the receiver

% Pr_rx = Pr_BD -