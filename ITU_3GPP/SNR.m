function [SINR_db] = SNR(f,d1,d2)

P_tx = 33; %dBm
% f = 200; %MHz

% lambda = (3*10^8)/(f*1e6);

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

% d3 = lambda:((10*lambda)-lambda)/10:(10*lambda); %multi-bounce link
% bd2bd_loss = D2D(f,d3); %BD to BD

% Losses in the link

fast_fading = 16; %dB
polarization_mismatch = 3; %dB
modulation_loss = 6; %dB
add_loss = fast_fading + polarization_mismatch + modulation_loss;

slow_fading_ITU = 6; %dB total 13, 7 dB is the SF in backscatter link (13-7 =6)
slow_fading_umi = 8.2;

%transmit antenna gain
G_t = 10;

%The value of P1 at BD1

UMi_loss = UMi(f,d1);
% P1_dbm = P_tx - (UMi_loss + add_loss + slow_fading_umi) + G_t;
% P1 = 10.^((P1_dbm-30)/10);

%% Loss due to ping-pong effect when the signal travels from Tx->BD1->BD2->Rx

% d2d_loss_pingpong = bd2bd_loss + modulation_loss; 
% rho = 10.^(-d2d_loss_pingpong/10);
% P_p = 0; % zeros(1,length(rho));
% for i=1:length(rho)
%     P_p(i) = P1*(rho(i)/(1-rho(i))); %accounts for all the ping-pong's
% end

%% Computation of P_rx (at the receiver)

d2d_loss = D2D(f,d2);
L_total_umi = d2d_loss + UMi_loss;

P_rx_dbm = P_tx - (L_total_umi + add_loss + slow_fading_umi) + G_t;
P_rx = 10.^((P_rx_dbm-30)/10);

%% SINR calculation

SINR_linear = P_rx/P_n;
% SINR_linear = ones(1,length(d1));
% for i = 1:length(d1)
%     SINR_linear(i) = P_rx/(P_n + P_p(i));
% end
SINR_db = 10*log10(SINR_linear);
end