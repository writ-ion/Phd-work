clc
clear all;
close all;

%% FSPL

f = [3.5*10^3 26*10^3];
d = 0:0.005:30;

L1 = 32.45 + 20*log10(f(1)) + 20*log10(d);
L2 = 32.45 + 20*log10(f(2)) + 20*log10(d);
figure(1)
plot(d,L1)
grid on
hold on
plot(d,L2)
xlabel('Distance (km)')
ylabel('Path loss (dB)')
legend ('3.5 GHz', '26 GHz')
hold off

%% Receiver sensitivity

% P_t_3.5 =  40 %in watts
% P_t_26 =  10 %in watts 5 W for small cells
B = [1e6 20e6 200e6];    %in Hz
NF = 8;       %noise figure in dB    
SNR = 4;     %in dB 
%B = B*10^6;   %bandwidth in hertz
T = 290;      %temperature in kelvin
k = 1.38064852*10^-23;     %Boltzmann constant
for i = B
    noise = 10*log10(k*T*B) + 30 + NF + SNR; %in dBm
    noise_linear = (10.^(noise/10)); %noise in linear scale
end
%% Radar equation

L_db = 0:1:20;
L = 10.^(L_db/10);       %Loss utilised from free space path loss
L = L.^(1/4);
Pt_dbm = 46;            %EIRP in dBm
Pt_lin = (1/1000)*10^(Pt_dbm/10);
Gr_dbi = 32; %receiver gain 32 dBi
Gt_dbi = 32; %transmitter gain (32 dBi)
Gt = 10^(Gt_dbi/10);                 
Gr = 10^(Gr_dbi/10);                 
c = 3*10^8;             %speed of light
% f = input('Enter the frequency (in GHz): ');
f = 3.5*10^9;           %operating frequency
lambda = c/f;
eta1 = 0.88*lambda^2
eta = [0.0004 eta1 0.01];  
% eta = [0.15 0.3 0.7];  
Pr_dbm = noise;
Pr = (1/1000)*10.^(Pr_dbm/10);

for i = 1:length(B)
    total_path_loss = Pt_dbm - (Pr_dbm) + (Gt_dbi + Gr_dbi)
end

%% Range calculation

r1 = ones(length(eta));
for i = 1:length(Pr_dbm)
    for j = 1:length(eta)
        r1(j,i) = ((Pt_lin*Gt*Gr*(lambda^2)*eta(j))/(((4*pi)^3)*Pr(i))).^(1/4);
    end
end
R = ones(length(L),length(eta)*3);
col_counter = 1;
row_counter=1;
for i = 1:length(eta)*3 %column counter
    for j = 1:length(L) %row counter
        R(j,i) = r1(i)/L(j); %distance in meters after additional loss
    end
end

R1=[R(:,1),R(:,2),R(:,3)];
R2=[R(:,4),R(:,5),R(:,6)];
R3=[R(:,7),R(:,8),R(:,9)];

figure(2)
%1 MHz
h1 = subplot(3,1,1); 
plot(R1(:,1),L_db,'-x','LineWidth',1.25);
hold on
plot(R1(:,2),L_db,'-.','LineWidth',1.25);
hold on
plot(R1(:,3),L_db,'-o','LineWidth',1.25);
% set(h1,{'Marker'},{'+';'s';'o';})
% legend('0.01','0.3','0.7')
title('1 MHz carrier bandwidth')
grid on

% 20 MHz
h2 = subplot(3,1,2);
plot(R2(:,1),L_db,'-x','LineWidth',1.25);
hold on
plot(R2(:,2),L_db,'-.','LineWidth',1.25);
hold on
plot(R2(:,3),L_db,'-o','LineWidth',1.25);
ylabel ('Additional loss (dB)');
title('20 MHz carrier bandwidth')
legend('0.0004','half-dipole','0.01')
grid on

% 200 MHz
h3 = subplot(3,1,3); 
plot(R3(:,1),L_db,'-x','LineWidth',1.25);
hold on
plot(R3(:,2),L_db,'-.','LineWidth',1.25);
hold on
plot(R3(:,3),L_db,'-o','LineWidth',1.25);
xlabel ('Distance (m)');
title('200 MHz carrier bandwidth')
grid on

% figure(3)
% h = plot(R,L_db,'LineWidth',1.25);
% grid on
% set(h,{'Marker'},{'+';'s';'o';'+';'s';'o';'+';'s';'o';})
% legend('1 MHz','20 MHz','200 MHz','1 MHz','20 MHz','200 MHz','1 MHz','20 MHz','200 MHz')
