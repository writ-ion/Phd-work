clear all
close all
clc

P_tx = 33; %dBm

f = 200; %MHz

d1_micro = 5:1.45:150; %15:1.35:150; %forward link, micro (150m maximum)
d1 = sqrt(d1_micro.^2 + (15-1.5)^2); %antenna height 1.5m

d2 = 0.6:0.6:60; % backscatter link

d3 = 0.1:0.1:10; % multi-bounce link
n = 1; %number of bounces, eg: 1 bounce is the two and fro signal 
...oscillation between 2 BDs

% Losses in the link
fast_fading = 16; %dB
polarization_mismatch = 3; %dB
modulation_loss = 6; %dB
add_loss = fast_fading + polarization_mismatch + modulation_loss;

slow_fading_ITU = 6; %dB total 13, 7 dB is the SF in backscatter link
...(13-7 =6)
% slow_fading_uma = 7.8;
slow_fading_umi = 8.2;

%transmit antenna gain
G_t = 10;

d2r_loss = D2D(f,d2); %BD to receiver
UMi_loss = UMi(f,d1,d1_micro); %forward link

for i = 1:length(d2r_loss)
    for j = 1:length(d2r_loss)
        %L_total_NLOS = d2d_loss + NLOS1(1:100);
        %L_total_uma = d2d_loss + UMa_loss(1:100);
        L_total_umi = d2r_loss + UMi_loss(1:100);
    end
end

%Ping-pong effect

d2d_loss = D2D(f,d3); %BD to BD

%% Case 1 : Signal from BD1 travels to BD2 and back

% one time modulation loss for the bounce from the active BD2 and two times 
% the path loss.
d2d_loss_firstbounce = 2 * d2d_loss + modulation_loss; %for one bounce

% after the signal returns to BD1 from BD2, additional modulation loss is
% added before the signal travels from BD1 to the RX.
d2d_loss_multibounce = (d2d_loss_firstbounce * n) + modulation_loss; %for n bounces

%% Case 2 : Signal from BD1 travels to BD2 and then further to RX

d2d_loss_singlebounce = n * (d2d_loss + modulation_loss);

%% Received signal level

rx_level_umi_case1 = P_tx + G_t - (L_total_umi + add_loss + slow_fading_umi)...
    - d2d_loss_multibounce(10);

rx_level_umi_case2 = P_tx + G_t - (L_total_umi + add_loss + slow_fading_umi)...
    - d2d_loss_singlebounce(10);

rx_level_umi_ref = P_tx + G_t - (L_total_umi + add_loss + slow_fading_umi);

%% Results and graphs

figure
mesh(d2,d1(1:100),rx_level_umi_case1)
grid minor
ylabel('Forward link [m]')
xlabel('Backscatter link [m]')
zlabel('Rx level [dBm]')
title('Case 1: Signal from Tx->BD1->BD2->BD1->Rx')

figure
mesh(d2,d1(1:100),rx_level_umi_case2)
grid minor
ylabel('Forward link [m]')
xlabel('Backscatter link [m]')
zlabel('Rx level [dBm]')
title('Case 2: Signal from Tx->BD1->BD2->Rx')

figure
mesh(d2,d1(1:100),rx_level_umi_ref)
grid minor
ylabel('Forward link [m]')
xlabel('Backscatter link [m]')
zlabel('Rx level [dBm]')
title('Reference: Signal from Tx->BD1->Rx')