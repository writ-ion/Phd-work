clc
clear all
close all

tic
f = 200;
lambda = (3*10^8)/(f*1e6);

%% location of the BD1, Tx and Rx

t_x = [0 0];
r_x = [100*lambda 0];

% Different locations for the BD1

theta = 135;
d = (10*lambda);

BD1 = [(r_x(1) + d*cos(theta)) (r_x(2) + d*sin(theta))];% 85.5 0;116 6.5;115 1.5];
s = size (BD1);
s = s(1);
% BD2 = [85.5 0];
% BD3 = [116 6.5];
% BD4 = [115 1.5];

% Location of the BD2

% X = 137;
% Y = 1.5;
% Z = [X; Y]';

X = 0:0.1:180;
Y = -50:0.1:50;
% X = 50:10:100;
% Y = -20:8:20;
[XX,YY] = ndgrid(X,Y);
Z = [XX(:),YY(:)];
% Z = [X; Y]';

%% Distance between Tx and BD1

d1 = ones(1,s);
i=1;
while i<= s
%for i = 1:length(BD1)
    d1(i) = distance(t_x(1,1),BD1(i,1),t_x(1,2),BD1(i,2));
    i = i+1;
%end
end
%% Distance between BD1 and Rx

d2 = ones(1,length(s));
for i = 1:length(s)
    d2(i) = distance(r_x(1),BD1(i,1),r_x(2),BD1(i,2));
end

%% Distance between Tx and BD2

d3 = ones(1,length(Z));
for i = 1:length(Z)
    d3(i) = distance(t_x(1),Z(i,1),t_x(2),Z(i,2));
end

%% Distance between BD2 and Rx

% d4 = ones(1,length(Z));
for i = 1:length(Z)
    d4(i) = distance(r_x(1),Z(i,1),r_x(2),Z(i,2));
end


%% Distance between BD1 and BD2

% bd2bd = ones(length(Z),length(BD1));
j=1;
while j <= length(s)
    for i = 1:length(Z)
        bd2bd(j,i) = distance(BD1(j,1),Z(i,1),BD1(j,2),Z(i,2));
    end
    j = j+1;
end

%% SINR and interference power calculation for all the BD1 points

SINR_CI = ones(1,length(Z));
SINR_SI = ones(1,length(Z));
SINR_TI = ones(1,length(Z));
Pi_cross = ones(1,length(Z));
Pi_self = ones(1,length(Z));
j=1;
while j <= length(s)
for i = 1:length(Z)
    [SINR_CI(j,i),SINR_SI(j,i),SINR_TI(j,i),Pi_cross(j,i),Pi_self(j,i)] = SINR_bd(f,d1(j),d2(j),d3(i),d4(i),bd2bd(j,i));
end
    j = j+1;
end
%
SINR_lim = -10; % in dB

SINR_cross = reshape(SINR_CI(1,:),length(X),length(Y));
SINR_cross(SINR_cross < SINR_lim) = NaN;
SINR_self = reshape(SINR_SI(1,:),length(X),length(Y));
SINR_self(SINR_self < SINR_lim) = NaN;
SINR_total = reshape(SINR_TI(1,:),length(X),length(Y));
SINR_total(SINR_total < SINR_lim) = NaN;
% SINR2 = reshape(SINR_BD(2,:),length(X),length(Y));
% SINR2(SINR2 < SINR_lim) = NaN;
% SINR3 = reshape(SINR_BD(3,:),length(X),length(Y));
% SINR3(SINR3 < SINR_lim) = NaN;
% SINR4 = reshape(SINR_BD(4,:),length(X),length(Y));
% SINR4(SINR4 < SINR_lim) = NaN;

%% Illustration of the environment

% figure
% scatter(t_x(1),t_x(2),'d');
% text(t_x(1),t_x(2),'Tx');
% hold on
% b1 = [50,50];
% b2 = [50,-50];
% plot(b1,b2,'-')
% % surf(XX,YY,SINR1);
% mesh(XX,YY,SINR_cross);
% colorbar
% plot(BD1(1),BD1(2),'r*','MarkerSize',20)
% scatter(BD1(1,1),BD1(1,2),'x');
% text(BD1(1,1),BD1(1,2),'BD1');
% scatter(r_x(1),r_x(2),'black','o','filled');
% text(r_x(1),r_x(2),'Rx');
% xlabel('link distance (m)')
% title('Cross Interference (limited to -10 dB)')
% axis([0 160 -50 50])
% grid on

% figure
% scatter(t_x(1),t_x(2),'d');
% text(t_x(1),t_x(2),'Tx');
% hold on
% b1 = [50,50];
% b2 = [50,-50];
% plot(b1,b2,'-')
% % surf(XX,YY,SINR1);
% mesh(XX,YY,SINR_self);
% colorbar
% plot(BD1(1),BD1(2),'r*','MarkerSize',20)
% scatter(BD1(1,1),BD1(1,2),'x');
% text(BD1(1,1),BD1(1,2),'BD1');
% scatter(r_x(1),r_x(2),'black','o','filled');
% text(r_x(1),r_x(2),'Rx');
% xlabel('link distance (m)')
% title('Self Interference (limited to -10 dB)')
% axis([0 160 -50 50])
% grid on

figure
% scatter(t_x(1),t_x(2),'d');
% text(t_x(1),t_x(2),'Tx');
hold on
% surf(XX,YY,SINR1);
mesh(XX,YY,SINR_total);
% colorbar
a=colorbar;
ylabel(a,'SINR (dB)','FontSize',16,'Rotation',270);
% plot(BD1(1),BD1(2),'r*','MarkerSize',2)
% 
% scatter(BD1(1,1),BD1(1,2),'x');
% text(BD1(1,1),BD1(1,2),'BD1');
% scatter(r_x(1),r_x(2),'black','o','filled');
% text(r_x(1),r_x(2),'Rx');
xlabel('X [m]')
ylabel('Y [m]')
% title('Cross Interference (limited to -10 dB)')
axis([0 160 -50 50])
grid on

%% CDF plot for SINR

% SINR_cdf = reshape(SINR,[],1);
% SINR_wb_cdf = reshape(SINR_wb,[],1);
% 
% figure
% h(1, 1) = cdfplot(SINR_cdf);
% hold on
% h(1, 2) = cdfplot(SINR_wb_cdf);
% set( h(:,1), 'LineStyle', ':', 'Color', 'r');
% set( h(:,2), 'LineStyle', '--', 'Color', 'b');
% xlabel('SINR (dB)')
% legend ('SINR - without bounce','SINR - with bounce')
% axis([18 19 0 1])

%% SINR plot for Interference power

CI_cdf = reshape(Pi_cross,[],1); %cross interference

Pi_self = SINR_self_BD2;
SI_cdf = reshape(Pi_self,[],1); %self interference

figure
h(1, 1) = cdfplot(CI_cdf);
hold on
h(1, 2) = cdfplot(SI_cdf);
set( h(:,1), 'LineStyle', '-', 'Color', 'r');
set( h(:,2), 'LineStyle', '--', 'Color', 'b');
xlabel('Interference Power [dBm]')
legend ('Cross Interference','Self Interference')
axis([-160 -100 0 1])

% figure
% cdfplot(CI_cdf)
% xlabel('Cross Interference Power (dBm)')
% 
% figure
% cdfplot(SI_cdf)
% xlabel('Self Interference Power (dBm)')

toc
%% function for calculating the distance between any two points

function [d] = distance(x1,x2,y1,y2)
    d = sqrt((x2-x1)^2 + (y2-y1)^2);
end