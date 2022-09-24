clc;
clear all;
close all;

%% Parameters

h_bs = 2;
h_ue = 1.5;

d = 1:0.1:11;
d_3d = sqrt(d.^2 + (h_bs-h_ue)^2);

f = 60; % 3GPP model
f1 = 60; % IEEE 802.11 ad model

Pt_dbm = 46; % EIRP in dBm
Pt_lin = (1/1000)*10^(Pt_dbm/10);

% Directivity of the antenna
U = Pt_lin/(4*pi); %radiation intensity


%% Indoor office environment (3GPP)

% LOS path loss
L_los = 32.4 + 17.3*log10(d_3d) + 20*log10(f) + 3; %shadow fading = 3 dB

% NLOS path loss
L_nlos = 17.3 + 38.3*log10(d_3d) + 24.9*log10(f) + 8.3; %shadow fading = 8.3 dB
% L_nlos = max(L_los,L_nlos);

% LOS probability for Mixed office
for i = 1:length(d_3d)
    if d_3d(i) < 1.2 || d_3d(i) == 1.2
        P_los(i) = 1;
    elseif d_3d(i) > 1.2 && d_3d(i) < 6.5
        P_los(i) = exp(-((d_3d(i)-1.2)/4.7));
    else
        P_los(i) = exp(-((d_3d(i)-6.5)/32.6))*0.32;
    end
    %Total loss after the probability is considered
        L_total_mixed(i) = P_los(i)*L_los(i) + (1-P_los(i))*L_nlos(i);
end

% LOS probability for Open office
for i = 1:length(d_3d)
    if d_3d(i) < 5 || d_3d(i) == 5
        P_los(i) = 1;
    elseif d_3d(i) > 5 && d_3d(i) < 49
        P_los(i) = exp(-((d_3d(i)-5)/70.8));
    else
        P_los(i) = exp(-((d_3d(i)-49)/211.7))*0.54;
    end
    %Total loss after the probability is considered
        L_total_open(i) = P_los(i)*L_los(i) + (1-P_los(i))*L_nlos(i);
end

figure
plot(d,L_los)
hold on
grid on
plot(d,L_nlos)
% plot(d_3d,L_total_mixed)
% plot(d_3d,L_total_open)
title('Indoor - Path loss vs distance')
% legend('LOS','NLOS','Total PL-mixed','Total PL-open')
xlabel('distance (m)');
ylabel('path loss (dB)');

%% 5GCM InH Office

% CI model with 1m reference distance
PL_los_5GCM = 32.4 + 17.3*log10(d_3d) + 20*log10(f) + 3.02; % shadow fading = 3.02 dB
plot(d,PL_los_5GCM)
% 5GCM single slope InH office

% CIF model
PL_CIF_singleSlope1 = 32.4 + 31.9*(1 + 0.06*((f-24.2)/24.2))*log10(d_3d) + 20*log10(f) + 8.29; % shadow fading = 8.29
plot(d,PL_CIF_singleSlope1,'s')

% ABG model
PL_ABG_singleSlope1 = 38.3*log10(d_3d) + 17.30 + 24.9*log10(f) + 8.03; %shadow fading = 8.03 dB
plot(d,PL_ABG_singleSlope1,'s')

% 5GCM dual slope InH Office

% CIF model
for i = 1:length(d)
    if d(i) < 7.8 || d(i) == 7.8
        PL_CIF_dualSlope1(i) = PL_los_5GCM(i) + 25.1*(1 + 0.12*((f-24.1)/24.1))*log10(d(i)) + 7.65; %shadow fading = 7.65
    else
        PL_CIF_dualSlope1(i) = PL_los_5GCM(i) + 25.1*(1 + 0.12*((f-24.1)/24.1))*log10(7.8) + 42.5*(1 + 0.04*((f-24.1)/24.1))*log10(d(i)/7.8) + 7.65; %shadow fading = 7.65
    end
end
plot(d,PL_CIF_dualSlope1,'s')

% ABG model
for i = 1:length(d)
    if d(i) < 6.9 || d(i) == 6.9
        PL_ABG_dualSlope1(i) = 17*log10(d(i)) + 33 + 24.9*log10(f) + 7.78; %shadow fading = 7.78
    else
        PL_ABG_dualSlope1(i) = 17*log10(6.9) + 33 + 24.9*log10(f) + 41.7*log10(d(i)/6.9) + 7.78; %shadow fading = 7.78
    end
end
plot(d,PL_ABG_dualSlope1,'s')

%% 5GCM InH Shopping mall

% CI model with 1m reference distance
PL_los_5GCM = 32.4 + 17.3*log10(d_3d) + 20*log10(f) + 1.18; %shadow fading = 1.18

% 5GCM single slope InH shopping mall

% CIF model
PL_CIF_singleSlope2 = 32.4 + 25.9*(1 + 0.01*((f-39.5)/39.5))*log10(d_3d) + 20*log10(f) + 7.4; % shadow fading = 7.4
plot(d,PL_CIF_singleSlope2,'o')

% ABG model
PL_ABG_singleSlope2 = 32.1*log10(d_3d) + 18.09 + 22.4*log10(f) + 6.97; %shadow fading = 6.97 dB
plot(d,PL_ABG_singleSlope2,'o')

% 5GCM dual slope InH Shopping mall

% CIF model
for i = 1:length(d)
    if d(i) < 110 || d(i) == 110
        PL_CIF_dualSlope2(i) = PL_los_5GCM(i) + 24.3*(1 - 0.01*((f-39.5)/39.5))*log10(d(i)) + 6.26; %shadow fading = 6.26
    else
        PL_CIF_dualSlope2(i) = PL_los_5GCM(i) + 24.3*(1 + 0.01*((f-39.5)/39.5))*log10(110) + 83.6*(1 + 0.39*((f-39.5)/39.5))*log10(d(i)/110) + 6.26; %shadow fading = 6.26
    end
end
plot(d,PL_CIF_dualSlope2,'o')

% ABG model
for i = 1:length(d)
    if d(i) < 147 || d(i) == 147
        PL_ABG_dualSlope2(i) = 29*log10(d(i)) + 22.17 + 22.4*log10(f) + 6.36; %shadow fading = 6.36
    else
        PL_ABG_dualSlope2(i) = 29*log10(147) + 22.17 + 22.4*log10(f) + 114.7*log10(d(i)/147) + 7.78; %shadow fading = 7.78
    end
end
plot(d,PL_ABG_dualSlope2,'o')

%% mmMAGIC InH Office

PL_los_mmMagic = 13.8*log10(d_3d) + 33.6 + 20.3*log10(f) + 1.18; %shadow fading = 1.18
plot(d,PL_los_mmMagic,'*')

PL_nlos_mmMagic = 36.9*log10(d_3d) + 15.2 + 26.8*log10(f) + 8.03; %shadow fading = 8.03
plot(d,PL_nlos_mmMagic,'*')

%% METIS InH Shopping Mall (frequency = 63 GHz)

PL_los_METIS = 68.8 + 18.4*log10(d) + 8.03; %shadow fading = 8.03, 1.5 < d < 13.4, h_bs = h_ ue = 2m
plot(d,PL_los_METIS,'d')

PL_nlos_METIS = 94.3 + 3.59*log10(d) + 8.03; %shadow fading = 8.03, 4 < d < 16.1, h_bs = h_ ue = 2m
plot(d,PL_nlos_METIS,'d')

%% IEEE 802.11 ad InH Office

PL_los_IEEE = 32.5 + 20*log10(f1) + 20*log10(d); %no shadow fading mentioned
plot(d,PL_los_IEEE,':')

PL_nlos_IEEE = 44.2 + 20*log10(f1) + 18*log10(d) + 1.5; %shadow fading = 1.5
plot(d,PL_nlos_IEEE,':')

%% figures

figure
plot(d,L_los)
hold on
plot(d,PL_los_5GCM)
plot(d,PL_los_mmMagic,'*')
plot(d,PL_los_METIS,'d')
plot(d,PL_los_IEEE,'-.')
grid on
xlabel('distance (m)');
ylabel('path loss (dB)');
title('LOS path loss')
legend('3GPP','5GCM','mmMAGIC','METIS','IEEE')

figure
plot(d,L_nlos)
hold on
plot(d,PL_CIF_singleSlope1,'-s')
plot(d,PL_ABG_singleSlope1,'-s')
plot(d,PL_CIF_dualSlope1,'s')
plot(d,PL_ABG_dualSlope1,'s')

plot(d,PL_CIF_singleSlope2,'-o')
plot(d,PL_ABG_singleSlope2,'-o')
plot(d,PL_CIF_dualSlope2,'o')
plot(d,PL_ABG_dualSlope2,'o')

plot(d,PL_nlos_mmMagic,'-d')
plot(d,PL_nlos_METIS,'.')
plot(d,PL_nlos_IEEE,'-.')
grid on
xlabel('distance (m)');
ylabel('path loss (dB)');
title('NLOS path loss')
legend('3GPP','5GCM-CIF single slope,office','5GCM-ABG single slope,office',...
    '5GCM-CIF dual slope,office','5GCM-ABG dual slope,office',...
    '5GCM-CIF single slope,shopping','5GCM-ABG single slope,shopping',...
    '5GCM-CIF dual slope,shopping','5GCM-ABG dual slope,shopping','mmMAGIC',...
    'METIS','IEEE')
hold off