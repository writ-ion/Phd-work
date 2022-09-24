clc;
close all;
clear all;

%% Parameters

L_db = 0:1:20;
L = 10.^(L_db/10);       %Loss utilised from free space path loss
L = L.^(1/4);
Pt_dbm = 36;            %EIRP in dBm
P_lin = (1/1000)*10^(Pt_dbm/10);
% Pt_dB = Pt_dbm - 30;    %dB to dBm
% Pt = 10^(Pt_dbm/10)    %transmit power in watts
Gr_dbi = 32; %receiver gain 32 dBi
Gt_dbi = 32; %transmitter gain (32 dBi)
Gt = 10^(Gt_dbi/10);                 
Gr = 10^(Gr_dbi/10);   
c = 3*10^8;             %speed of light
f = 26*10^6;           %operating frequency
lambda = c/f;
eta1 = 0.88*lambda^2;
eta = [0.001 eta1 0.01 0.15 0.3 0.7];              %radar non-fluctuating cross section (in meters)(.88*lambda^2)half dipole antenna
% Rt = 100;             %Transmitter to target = 100 meters
% Rr = 100;             %Receiver to target = 100 meters (monostatic radar)
Pr_dbm = -93.97;
Pr = (1/1000)*10^(Pr_dbm/10);
%R = 1560;                %Range from Tx to sensor (monostatic radar)

%% Range calculation

r1 = ((P_lin*Gt*Gr*(lambda^2)*eta)/(((4*pi)^3)*Pr)).^(1/4);
% R1 = r1/1000;
% r2 = ((P_lin*Gt*Gr*(lambda^2)*eta)/(((4*pi)^3)*Pr*L0)).^(1/4);
% R2 = r2/1000;
R=ones(length(L),length(eta));
col_counter = 1;
row_counter=1;
while col_counter<6
    for range = r1
        row_counter=1;
        for i=L
            R(row_counter,col_counter) = range/(i*1000);
            row_counter=row_counter+1;
        end
        col_counter=col_counter+1;
    end
end
h = plot(R,L_db,'LineWidth',1.25);
set(h,{'Marker'},{'+';'s';'o';'*';'d';'pentagram'})
% plot(R,L_db)
grid on
% axis([0 11 0 20])
legend('0.001','0.01','0.15','0.3','0.7','7.92')
ylabel('Additional loss (dB)')
xlabel('Distance (km)')

% R = [R1;R2];
% bar(R')
% grid on
% legend('show')
