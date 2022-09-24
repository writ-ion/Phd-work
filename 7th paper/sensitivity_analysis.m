clc
clear all
close all

%% Receiver sensitivity

k = 1.38*10^-23;
T = 290;
B = 12*15; % 12RB*15khz (LTE)
B = B*10^3;

interference_ideal = 20;     % ideal cell
interference_loaded = -10;   % loaded cell
NF = 10;

thermal_noise = 10*log10(k*T*B/0.001) + NF;         %minimum signal level

Pr_ideal = thermal_noise + interference_ideal;      %noise floor for unloaded cell (interference)
Pr_loaded = thermal_noise + interference_loaded;    %noise floor for loaded cell

%% Parameters

Pt = 46;                                    % transmit power in dBm
P_lin = (1/1000)*10^(Pt/10);                % transmit power in linear scale

Gt_dbi = 18;                                %transmitter gain (18 dB)
Gt = 10^(Gt_dbi/10);   
Gr_dbi = 0;                                 %receiver gain
Gr = 10^(Gr_dbi/10);
              
Pr_ideal_lin = (1/1000)*10^(Pr_ideal/10);   %noise floor of unloaded cell in linear scale
Pr_loaded_lin = (1/1000)*10^(Pr_loaded/10); %noise floor of loaded cell in linear scale

d = 0.5:1:50;                               % distance in kilometers
f = 700;                                    % frequency in megahertz
lambda = (3*10^8)/(f*10^6);                 % wavelength in meters

eta = 0.88*lambda^2;                        % radar cross section from Barton's book

L_db = 5;                                   % additional loss (in dB)
L_lin = 10^(L_db/10);                       % additional loss (in linear sccale)

%% Free space path loss

FSPL = 32.4 + 20*log10(d) + 20*log10(f);
plot (d,FSPL)
grid on

available_loss_ideal = Pt - Pr_ideal;
available_loss_loaded = Pt - Pr_loaded;

%% Radar equation

R_ideal = ((P_lin*Gt*Gr*(lambda^2)*eta)/(((4*pi)^3)*Pr_ideal_lin*L_lin)).^(1/4)

R_loaded = ((P_lin*Gt*Gr*(lambda^2)*eta)/(((4*pi)^3)*Pr_loaded_lin*L_lin)).^(1/4)