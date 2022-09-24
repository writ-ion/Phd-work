clc;
clear all;
close all;

%% Simulation parameters
load 'rss.mat';
f=100; %fm frequency MHz
Tx_pow=60000; %60kW
P_dbm=10*log10(Tx_pow)+30;

Tx1=30; %receiver height
Tx2=80;
Rx=1; %sensor height meters 

diffraction_loss=20; %dB
sensor_loss=10; %dB
coupling_loss=40;

min_reception_level=-123.94; %dBm

%% path loss from sensor to Rx antenna
L1=ones(1,length(rss));
L2=ones(1,length(rss));
L3=ones(1,length(rss));
L4=ones(1,length(rss));

for i=1:length(rss)
    L1(1,i)=-min_reception_level+(rss(1,i)-sensor_loss-diffraction_loss); %Okumura Hata
    L2(1,i)=-min_reception_level+(rss(2,i)-sensor_loss-diffraction_loss); %FSPL
    L3(1,i)=-min_reception_level+(rss(1,i)-coupling_loss);
    L4(1,i)=-min_reception_level+(rss(2,i)-coupling_loss);
end

%link budget formula
%received power(dBm)=Transmitted power(dBm)+gains(dB)-losses(dB)

%% Cost 231 Hata model parameters

A=69.55;
B=26.16;
C=39.5; %propagation slope reduced from 44.9 to ensure > 2.5
Cm=-10; %area correction factor open rural area

%large city
a_hms1=3.2*(log10(11.75*Rx))^2-4.97; 

%% distance
% Okumura Hata

d1=ones(1,length(L1));
d2=ones(1,length(L1));
d3=ones(1,length(L1));
d4=ones(1,length(L1));
for i=1:length(L1)
    d1(1,i) = 10^((L1(1,i)-A-B*log10(f)+13.82*log10(Tx1)+a_hms1-Cm)/(C-6.55*log10(Tx1))); %sensor loss
    d2(1,i) = 10^((L1(1,i)-A-B*log10(f)+13.82*log10(Tx2)+a_hms1-Cm)/(C-6.55*log10(Tx2)));
    d3(1,i) = 10^((L3(1,i)-A-B*log10(f)+13.82*log10(Tx1)+a_hms1-Cm)/(C-6.55*log10(Tx1))); %coupling loss considered
    d4(1,i) = 10^((L3(1,i)-A-B*log10(f)+13.82*log10(Tx2)+a_hms1-Cm)/(C-6.55*log10(Tx2)));
end

% Free space path loss

d_fspl=ones(1,length(L2));
d_fspl1=ones(1,length(L2));
for j=1:length(L2)
    d_fspl(1,j) = 10^((L2(1,j)-32.45-20*log10(f))/20);  %sensor loss
    d_fspl1(1,j) = 10^((L4(1,j)-32.45-20*log10(f))/20); %coupling loss
end