clc
clear all
close all

P_tx_macro = 43; %dBm
P_tx_micro = 33;

f = 200; %MHz

d1_micro = 5:1.45:150; %15:1.35:150; %forward link, micro (150 maximum) 5
d1_micro_actual = sqrt(d1_micro.^2 + (15-1.5)^2);

d1_macro = 5:2.45:250; %30:2.2:250; % forward link, macro (250 maximum) 25
d1_macro_actual = sqrt(d1_macro.^2 + (30-1.5)^2);

d2 = 0.6:0.6:60; % backscatter link

% Losses in the link
fast_fading = 16; %dB
polarization_mismatch = 3; %dB
modulation_loss = 6; %dB
add_loss = fast_fading + polarization_mismatch + modulation_loss;

slow_fading_ITU = 6; %dB total 13, 7 dB is the SF in backscatter link (13-7 =6)
slow_fading_uma = 7.8;
slow_fading_umi = 8.2;

%transmit antenna gain
G_t = 10;

% RX level of LoRa
rx_lora = -141; %dBm

d2d_loss = D2D(f,d2);
NLOS1 = NLOS1(f,d1_macro_actual);
UMa_loss = UMa(f,d1_macro_actual);
UMi_loss = UMi(f,d1_micro_actual);

for i = 1:length(d2d_loss)
    for j = 1:length(d2d_loss)
        L_total_NLOS = d2d_loss + NLOS1(1:100);
        L_total_uma = d2d_loss + UMa_loss(1:100);
        L_total_umi = d2d_loss + UMi_loss(1:100);
    end
end

rx_level_NLOS = P_tx_macro - (L_total_NLOS + add_loss + slow_fading_ITU) + G_t;

rx_level_uma = P_tx_macro - (L_total_uma + add_loss + slow_fading_uma) + G_t;
rx_level_umi = P_tx_micro - (L_total_umi + add_loss + slow_fading_umi) + G_t;

%Assigning NULL for values exceding the min RX level for LoRa

% rx_level_NLOS(rx_level_NLOS < rx_lora) = NaN;
% 
% rx_level_uma(rx_level_uma < rx_lora) = NaN;
%rx_level_umi(rx_level_umi < rx_lora) = NaN;

%Plotting of the graphs
d = d2 + d1_macro_actual(1:100); %total distance for the 2D graphs

% figure
% mesh(d2,d1_macro_actual(1:100),rx_level_NLOS)
% grid minor
% ylabel('Forward link (m)')
% xlabel('Backscatter link (m)')
% zlabel('Rx level (dBm)')
% title('ITU, rooftop-to-street')

% figure
% plot(d,rx_level_NLOS)
% grid minor
% xlabel('Forward + Backscatter link (m)')
% ylabel('Rx level (dBm)')
% title('ITU, rooftop-to-street')
% axis([min(d) max(d) min(min(rx_level_NLOS)) max(max(rx_level_NLOS))])

% figure
% mesh(d2,d1_macro_actual(1:100),rx_level_uma)
% grid minor
% ylabel('Forward link (m)')
% xlabel('Backscatter link (m)')
% zlabel('Rx level (dBm)')
% title('3GPP, Urban macrocell')

% figure
% plot(d,rx_level_uma)
% grid minor
% xlabel('Forward + Backscatter link (m)')
% ylabel('Rx level (dBm)')
% title('3GPP, Urban macrocell')
% axis([min(d) max(d) min(min(rx_level_uma)) max(max(rx_level_uma))])

figure
mesh(d2,d1_micro_actual(1:100),rx_level_umi)
grid minor
ylabel('Forward link [m]')
xlabel('Backscatter link [m]')
zlabel('Rx level [dBm]')
title('3GPP, Urban microcell')

% figure
% plot(d,rx_level_umi)
% grid minor
% xlabel('Forward + Backscatter link (m)')
% ylabel('Rx level (dBm)')
% title('3GPP, Urban microcell')
% axis([min(d) max(d) min(min(rx_level_umi)) max(max(rx_level_umi))])
