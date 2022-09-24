%% noise calculations

P_dbm = 46; %dBm
% add_loss_suburban = 5; % in dB
% add_loss_urban = 15;
B = 12*15;    %in kHz
NF = 10;       %noise figure in dB    
SNR = 2;     %in dB 
B= B*10^3;   %bandwidth in hertz
T = 290;      %temperature in kelvin
k = 1.38064852*10^-23;     %boltzmann constant
% noise1 = 10*log10((k*T*B)/0.001) + NF + SNR
noise = 10*log10(k*T*B) + 30 + NF + SNR; %in dBm
noise_linear = (10^(noise/10)); %noise in linear scale
path_loss = P_dbm - noise
% total_loss_suburban = path_loss - add_loss_suburban
% total_loss_urban = path_loss - add_loss_urban

% f = 3500;
% d = 5;
% L = 34.0576 + 31*log10(f) + 23*log10(d)