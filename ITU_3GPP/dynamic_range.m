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

% The following calculates the RX level utilising the different parameters
% for the direct link from the TX to the RX, the distance is d1+d2 in this
% case

for i = 1:length(d2)
    for j = 1:length(d2)
        d_macro(i,j) = d1_macro_actual(i) + d2(j);
        d_micro(i,j) = d1_micro_actual(i) + d2(j);
    end
end

for i = 1:length(d2)
    for j = 1:length(d2)
%         d_macro(i,j) = d1_macro_actual(i) + d2(j);
        d_micro_new(i,j) = d1_micro(i) + d2(j);
    end
end

% d_macro = d1_macro_actual + d2;
% d_micro = d1_micro_actual + d2;
for i = 1:length(d2)
    for j = 1:length(d2)
        L_NLOS1(i,j) = NLOS1(f,d_macro(i,j));
        L_uma_loss(i,j) = UMa(f,d_macro(i,j));
        L_umi_loss(i,j) = UMi(f,d_micro(i,j),d_micro_new(i,j));
    end
end

rx_level_NLOS_direct = P_tx_macro - (L_NLOS1 + add_loss + slow_fading_ITU) + G_t;

rx_level_uma_direct = P_tx_macro - (L_uma_loss + add_loss + slow_fading_uma) + G_t;
rx_level_umi_direct = P_tx_micro - (L_umi_loss + add_loss + slow_fading_umi) + G_t;

% The following computes the RX level utilising the different parameters
% for the backscatter communication, where d1 and d2 are separate links,
% between the Tx and sensor, sensor and Rx.

d2d_loss = D2D(f,d2);
NLOS1 = NLOS1(f,d1_macro_actual);
UMa_loss = UMa(f,d1_macro_actual);
UMi_loss = UMi(f,d1_micro_actual,d1_micro);

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

% figure
% mesh(d2,d1_micro_actual(1:100),rx_level_NLOS_direct)
% hold on 
% mesh(d2,d1_micro_actual(1:100),rx_level_NLOS)
% grid minor
% ylabel('Forward link (m)')
% xlabel('Backscatter link (m)')
% zlabel('Rx level (dBm)')
% title('Dynamic range comparison - ITU, rooftop-to-street')
% legend('Direct','Backscatter')
% 
% NLOS_rx_level = abs(rx_level_NLOS_direct - rx_level_NLOS);
% figure
% mesh(d2,d1_micro_actual(1:100),NLOS_rx_level)
% grid minor
% ylabel('Forward link (m)')
% xlabel('Backscatter link (m)')
% zlabel('Rx level (dBm)')
% title('Dynamic range - ITU, rooftop-to-street')
% % legend('Direct','Backscatter')
% 
% figure
% mesh(d2,d1_micro_actual(1:100),rx_level_uma_direct)
% hold on 
% mesh(d2,d1_micro_actual(1:100),rx_level_uma)
% grid minor
% ylabel('Forward link (m)')
% xlabel('Backscatter link (m)')
% zlabel('Rx level (dBm)')
% title('Dynamic range comparison - 3GPP, Urban macrocell')
% legend('Direct','Backscatter')
% 
% uma_rx_level = abs(rx_level_uma_direct - rx_level_uma);
% figure
% mesh(d2,d1_micro_actual(1:100),uma_rx_level)
% grid minor
% ylabel('Forward link (m)')
% xlabel('Backscatter link (m)')
% zlabel('Rx level (dBm)')
% title('Dynamic range - 3GPP, Urban macrocell')
% % legend('Direct','Backscatter')
% 
% figure
% mesh(d2,d1_micro_actual(1:100),rx_level_umi_direct)
% hold on 
% mesh(d2,d1_micro_actual(1:100),rx_level_umi)
% grid minor
% ylabel('Forward link (m)')
% xlabel('Backscatter link (m)')
% zlabel('Rx level (dBm)')
% title('Dynamic range comparison - 3GPP, Urban microcell')
% legend('Direct','Backscatter')

umi_rx_level = abs(rx_level_umi_direct - rx_level_umi);
% umi_rx_level(umi_rx_level > 30) = NaN;
figure
mesh(d2,d1_micro_actual(1:100),umi_rx_level)
grid minor
ylabel('Forward link [m]')
xlabel('Backscatter link [m]')
zlabel('Dynamic range [dB]')
title('Dynamic range - Urban microcell')
% legend('Direct','Backscatter')