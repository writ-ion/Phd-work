clc
clear all
close all
tic

f = 200;
lambda = (3*10^8)/(f*1e6);

%% location of the BD1, Tx and Rx

t_x = [0 0];
r_x = [100*lambda 0];


% theta = 45;
% x = 20*lambda*cos(theta);
% y = 20*lambda*sin(theta);
X = 0:0.1:180;
Y = -50:0.1:50;
[XX,YY] = ndgrid(X,Y);
Z = [XX(:),YY(:)];
%Z = [X; Y]';

%% Distance between Tx and BD1 
d1 = ones(1,length(Z));
for i = 1:length(Z)
    d1(i) = distance(t_x(1),Z(i,1),t_x(2),Z(i,2));
end

%% Distance between BD1 and Rx
d2 = ones(1,length(Z));
for i = 1:length(Z)
    d2(i) = distance(r_x(1),Z(i,1),r_x(2),Z(i,2));
end

%% SINR calculation for all the BD1 points

SNR_BD = ones(1,length(Z));
for i = 1:length(Z)
    SNR_BD(i) = SNR(f,d1(i),d2(i));
end

SNR = reshape(SNR_BD,length(X),length(Y));
SNR_lim = -10;
SNR(SNR < SNR_lim) = NaN;

%% Illustration of the environment when only 1 backscatter device is present

figure
scatter(t_x(1),t_x(2),'d');
text(t_x(1),t_x(2),'Tx');
hold on
% scatter(r_x(1),r_x(2),'o');
% text(r_x(1),r_x(2),'Rx');
%plot(Z(:,1),Z(:,2),'rx')
% b1 = [50,50];
% b2 = [50,-50];
% plot(b1,b2,'-')
% surf(XX,YY,SINR);
mesh(XX,YY,SNR);
colorbar
xlabel('relative link distance (m)')
title('SNR (limited to -10) when only one BD is present in the network')
axis([0 160 -50 50])
grid on
toc
%% function for calculating the distance between any two points

function [d] = distance(x1,x2,y1,y2)
    d = sqrt((x2-x1)^2 + (y2-y1)^2);
end