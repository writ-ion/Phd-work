clc
clear all
close all

%% Parameters

h_bs = 2; % height of base station
h_ue = 1.5; % height of user equipment

d1 = 0:0.1:50; % linear distance between Tx and RIS (in meters)
d2 = 2; % linear distance between RIS and Rx (in meters)

d = 2:0.1:52; %total distance between the TX and RX

% d_3d = sqrt(d1.^2 + (h_bs-h_ue)^2);

f = 60; % frequency (in GHz)
l = (3*10^8)/(f*10^9); % wavelength (in meters)

Pt_dbm = 33; % Transmit power in dBm
Pt = (1/1000)*10^(Pt_dbm/10);

% shadowing = 1.5;
% small_scale_fading = 6; % in dB

degree = 0;
dir = degree*pi/180;
syms theta;

%Calculation of the directive gain

P_rad = double(int(sin(theta)*cos(theta),theta,[0 pi/2])); % cos(theta)^3,int(sin(theta)*cos(theta))

D = (4*pi)/(2*pi)*P_rad; %antenna directivity
E = 1; %antenna efficiency + other losses including polarization mismatch and modulation losses
G = E*D;

%Gt = 27; % Gain of the transmit antenna
Gt_db = 27; % Gain of the transmit antenna (in dBi)
% Gt = 10^(Gt_db/10);
Gt = 1;

Gr = 1; % Gain of the receiver antenna

% G = 1; % Gain of the unit cell
A = 1; % reflection gain (considered unity here)
dx = 0.01; % size of unit cell in x-axis (in meters)
dy = 0.01; % size of unit cell in y-axis (in meters)
M = 8; % number of rows for the regularly arranged RIS
N = 8; % number of columns for the regularly arranged RIS

bound = (2*M*N*dx*dy)/l;
%% Received signal power in RIS-assisted wireless communications

%far field beamforming case 
% (The near field lower bound in the paper is considered to be 5*lambda
% and any distance greater than this specified distance is considered...
% to be in the far field of the Tx/Rx.

Pr = (Pt*Gt*Gr*G*(M^2)*(N^2)*dx*dy*(l^2)*(A^2))./(64*(pi^3)*(d1.^2)*d2^2); % Received power

Pr_comb = (Pt*Gt*Gr*G*(M^2)*(N^2)*dx*dy*(l^2)*(A^2))./(64*(pi^3)*(d.^2)); % Received power of the combined distance


PL = (64*(pi^3)*(d1.^2)*d2^2)./(Gt*Gr*G*(M^2)*(N^2)*dx*dy*(l^2)*(A^2)); % Path loss

PL_comb = (64*(pi^3)*(d.^2))./(Gt*Gr*G*(M^2)*(N^2)*dx*dy*(l^2)*(A^2)); % Path loss of the combined distance

% here the N and M terms are even numbers and represent the number of
% RIS elements and the dimensions of each RIS element is given by dx and
% dy, respectively. d1 & d2 represent the distance of thr RIS from the TX
% and the RX respectively and A represents the unit gain of the individual
% RIS element.

% figure
% plot(d1,(1./PL))
% grid on
% xlabel('distance (m)')
% ylabel('path loss (dB)')

figure
plot(d1,(10*log10(PL)))
hold on
plot(d1,(10*log10(PL_comb)))
grid on
xlabel('distance (m)')
ylabel('path loss (dB)')
legend('d1+d2','d')

figure
plot(d1,10*log10(Pr))
hold on
plot(d1,(10*log10(Pr_comb)))
grid on
xlabel('distance (m)')
ylabel('received power')
legend('d1+d2','d')