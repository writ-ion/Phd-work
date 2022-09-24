clc;
clear all;
close all;

R1 = 30; % distance from Tx to sensor in m (I have fixed this)
% R2=1:50; %distance from Tx to sensor in km
R = 0:100; % total link distance
%R=R1+R2; %distance from Tx to Rx

y=30; %additional path loss
Y=10^(y/10); %in linear scale

c=3e8; %in Hertz
f=100e6; %in Hertz
lambda=c/f; %wavelength

sigma = 1;%.88*lambda^2; %0.1:0.1:1;
DR = 5:0.5:55; %dynamic range in dB scale
dr = 10.^(DR/10);

%% Riku Jäntti's Calculations

temp = sqrt(sigma.*dr/(4*pi));
R2=(temp./R1)'.*R;

mesh(DR,R,R2)

% % R2 = (temp.*R')./ R1;
% R2 = R1'*(sqrt(sigma/(4*pi)*DR))./R;
% mesh(R2)    
xlabel('Dynamic range (dB)')
ylabel('Link distance (km)')
zlabel('Path after backscatterer(km)')
title('Professor Jäntti calculations (R1=30km)')

%% when the backscatter link is a simple sum of R1 & R2 = (R1+R2)

X=R'.*sqrt(dr)-R1-Y; %calculating R2, here represented by X
figure;
mesh(DR,R,X)
axis([0 55 0 100 0 600])
xlabel('Dynamic range (dB)')
ylabel('Link distance (km)')
zlabel('Path after backscatterer(km)')
title('Backscatter link is R1(=30km)+R2')

%% New calculations with path loss exponent of 2.4

BC= (((R.^2)'.*dr)*1000).^(1/2.4) - R1; % when path loss is +30 dB
at_least = 0;   %for example km
at_most = 300;     %for example km
new_matrix = min(max(BC, at_least), at_most);
figure;
mesh(DR,R,new_matrix)
xlabel('Dynamic range (dB)')
ylabel('Link distance (km)')
zlabel('Path after backscatterer(km)')
title('Backscatter link is R1(=30km)+R2+30 dB additional')

%% New calculations with path loss exponent of 2.4

% BC = (((R.^2)'.*dr)/1000).^(1/2.4) - R1; % when path loss is -30 dB
% at_least = 0;   %for example km
% at_most = 300;     %for example km
% new_matrix = min(max(BC, at_least), at_most);
% figure;
% mesh(DR,R,new_matrix)
% xlabel('Dynamic range (dB)')
% ylabel('Link distance (km)')
% zlabel('Path after backscatterer(km)')
% title('Backscatter link is R1(=30km)+R2-30 dB additional')

%% when the backscatter link is a combination of two separate paths, R1 & R2

p = sqrt(dr./R1^2-(dr'.*(R.^2))) ;%result is complex valued
q = R.*R1;
pq = q.*p;
% mesh(DR,R,pq)