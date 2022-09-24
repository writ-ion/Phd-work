function [SINR_db_cross,SINR_db_self,SINR_db_total,Pi_cross_log,Pi_self_log] = SINR_bd(f,d1,d2,d3,d4,bd2bd)
%
P_tx = 33; %dBm
% f = 60; %MHz

lambda = (3*10^8)/(f*1e6);

% Thermal noise calculation

B = 10;    %in kHz
B = B*1e3;   %bandwidth in hertz
T = 290;      %temperature in kelvin
k = 1.38064852 * 10^-23;     %boltzmann constant

noise = 10*log10((k*T*B)/0.001); %in dBm
P_n = (10^((noise-30)/10)); %noise in linear scale

%Distance description of the direct, forward and backscatter and multibounce links

% d1_micro = 20*lambda; %forward link
% d1 = d1_micro/sqrt(2);
% 
% d2 = 10*lambda; % backscatter link

%d3 = lambda:((10*lambda)-lambda)/10:(10*lambda); %multi-bounce link


% Losses in the link

fast_fading = 16; %dB
polarization_mismatch = 3; %dB
modulation_loss = 6; %dB
add_loss = fast_fading + polarization_mismatch + modulation_loss;

slow_fading_ITU = 6; %dB total 13, 7 dB is the SF in backscatter link (13-7 =6)
slow_fading_umi = 8.2;

%transmit antenna gain
G_t = 10;

%% Calculation of path loss for each individual link (d1,d2,d3,d4 and bd2bd)

d1_loss = UMi(f,d1);
d2_loss = D2D(f,d2);
d3_loss = UMi(f,d3);
d4_loss = D2D(f,d4);
bd2bd_loss = D2D(f,bd2bd) + modulation_loss;

%The value of P1 at BD1
% P1_dbm = P_tx - (UMi_loss + add_loss + slow_fading_umi) + G_t;
% P1 = 10^((P1_dbm-30)/10);

%% Computation of P_rx (at the receiver for BD1, path is Tx->BD1->Rx)

L_total_d1_d2 = d2_loss + d1_loss; %path loss for d1 and d2

P_rx_dbm = P_tx - (L_total_d1_d2 + add_loss + slow_fading_umi) + G_t;
P_rx = 10.^((P_rx_dbm-30)/10); % signal power

%% Loss due to ping-pong effect when the signal travels from Tx->BD1->BD2->Rx

% bd2bd_loss = D2D(f,bd2bd); %BD to BD
% bd2bd_loss = bd2bd_loss + modulation_loss; 
% rho = 10.^(-d2d_loss_pingpong/10);
% P_p = zeros(1,length(rho));
% for i=1:length(rho)
%     P_p(i) = P1*(rho(i)/(1-rho(i))); %accounts for all the ping-pong's
% end

%% Computation of Self Interference at the receiver (Tx->BD1->BD2->Rx)

L_total_interference1 = d1_loss + bd2bd_loss + d4_loss;

P_rx_interference1_dbm = P_tx - (L_total_interference1 + add_loss + slow_fading_umi) + G_t;
P_i1 = 10.^((P_rx_interference1_dbm-30)/10);

%% Computation of Self Interference at the receiver, 2 bounces (Tx->BD1->BD2->BD1->Rx)

L_total_interference12 = d1_loss + 2*bd2bd_loss + d2_loss;

P_rx_interference1_dbm = P_tx - (L_total_interference12 + add_loss + slow_fading_umi) + G_t;
P_i12 = 10.^((P_rx_interference1_dbm-30)/10);

%% Computation of Cross Interference at the receiver (Tx->BD2->BD1->Rx)

L_total_interference2 = d3_loss + bd2bd_loss + d2_loss;

P_rx_interference2_dbm = P_tx - (L_total_interference2 + add_loss + slow_fading_umi) + G_t;
P_i2 = 10.^((P_rx_interference2_dbm-30)/10);

%% Computation of Cross Interference at the receiver (Tx->BD2->Rx)

L_total_interference3 = d3_loss + d4_loss;

% P_rx_dbm = P_tx - (L_total_d1_d2 + add_loss + slow_fading_umi) + G_t;

P_rx_interference3_dbm = P_tx - (L_total_interference3 + add_loss + slow_fading_umi) + G_t;
P_i3 = 10.^((P_rx_interference3_dbm-30)/10);

%% Interference power calculation

Pi_cross = ones(1,length(d3));
Pi_self = ones(1,length(d3));
for i = 1:length(d3)
    Pi_cross(i) = P_i2(i) + P_i3(i);
    Pi_self(i) = P_i1(i);
end

Pi_cross_log = 10*log10(Pi_cross)+30;
Pi_self_log = 10*log10(Pi_self)+30;

%% SINR calculation

SINR_linear_cross = ones(1,length(d3));
SINR_linear_self = ones(1,length(d3));
SINR_linear_total = ones(1,length(d3));
for i = 1:length(d3)
    SINR_linear_cross(i) = P_rx./(P_n + (P_i2(i)));
    SINR_linear_self(i) = P_rx./(P_n + P_i1(i));
    SINR_linear_total(i) = P_rx./(P_n + (P_i2(i) + P_i3(i)));
end
SINR_db_cross = 10*log10(SINR_linear_cross);
SINR_db_self = 10*log10(SINR_linear_self);
SINR_db_total = 10*log10(SINR_linear_total);
%imagesc(SINR_db)
end