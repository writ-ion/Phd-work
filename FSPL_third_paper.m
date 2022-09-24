clc;
clear all;
close all;

%% Simulation parameters
f=100; %fm frequency MHz
P_dbm = 77.78;
d = 0:0.5:100;

%link budget formula 
%received power(dBm)=Transmitted power(dBm)+gains(dB)-losses(dB)

%%FSPL
L=32.45+20*log10(d)+20*log10(f);

% Pr_dbm = P_dbm + 32 - L;
% plot(d,Pr_dbm,'-*');
% grid on;
% xlabel('Distance (km)')
% ylabel('Received signal strength (dBm)')
% title('Signal strength vs distance')

% Path Loss

figure;
plot(d,L,'-*');
grid on;
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
title('Path loss vs distance');

%% Receiver sensitivity

B = 1e3;    %in Hz
NF = 10;       %noise figure in dB    
SNR = 10;     %in dB 
%B = B*10^6;   %bandwidth in hertz
T = 290;      %temperature in kelvin
k = 1.38064852*10^-23;     %Boltzmann constant
noise = 10*log10(k*T*B) + 30 + NF + SNR; %in dBm
noise_linear = (10.^(noise/10)); %noise in linear scale

%% Parameters

L_db = 10:1:30;
L = 10.^(L_db/10);       %Loss utilised from free space path loss
L = L.^(1/4);

Pt_dbm = 77.78;            %EIRP in dBm
P_lin = (1/1000)*10^(Pt_dbm/10);
% Pt_dB = Pt_dbm - 30;    %dB to dBm
% Pt = 10^(Pt_dbm/10)    %transmit power in watts

Gt_dbi = 0;                 %transmitter gain (16 dB)
Gr_dbi = 0;                 %receiver gain
Gt = 10^(Gt_dbi/10);                 
Gr = 10^(Gr_dbi/10);
% Gt = 1;
% Gr = 1;

c = 3*10^8;             %speed of light
f = 100*10^6;           %operating frequency
lambda = c/f;

sigma = 0.88*(lambda^2)

eta = [0.001 0.01 0.16 0.3 0.7 sigma];  %radar non-fluctuating cross section (in meters)

Pr_dbm = noise;
Pr = (1/1000)*10^(Pr_dbm/10);

% cable_loss = 2;
% total_path_loss = Pt_dbm - (Pr_dbm + cable_loss) + (Gt_dbi + Gr_dbi)
% 
% d = 10^((total_path_loss - 32.45 + 20*log10(f))/20)

%% Range calculation

r1 = ((P_lin*Gt*Gr*(lambda^2)*eta)/(((4*pi)^3)*Pr)).^(1/4);
% R1 = r1/1000;
% r2 = ((P_lin*Gt*Gr*(lambda^2)*eta)/(((4*pi)^3)*Pr*L0)).^(1/4);
% R2 = r2/1000;
R=ones(length(L),length(eta));
col_counter = 1;
row_counter=1;
while col_counter<=5
    for range = r1
        row_counter=1;
        for i=L
            R(row_counter,col_counter) = range/(i*1000);
            row_counter=row_counter+1;
        end
        col_counter=col_counter+1;
    end
end
figure
h = plot(R,L_db,'LineWidth',1.25);
set(h,{'Marker'},{'+';'s';'o';'*';'d';'<'})
% plot(R,L_db)
grid on
% grid minor
% axis([0 11 0 20])
legend('0.001','0.01','0.16','0.3','0.7','7.92')
ylabel('Additional loss (dB)')
xlabel('Distance (km)')