clc
clear all
close all

f = 200;
lambda = (3*10^8)/(f*1e6);

%% location of the BD1, Tx and Rx

t_x = [0 0];
r_x = [100 0];

% theta = 45;
% x = 20*lambda*cos(theta);
% y = 20*lambda*sin(theta);
X = 50:10:100;
Y = -20:8:20;
[XX,YY] = ndgrid(X,Y);
Z = [XX(:),YY(:)];
%Z = [X; Y]';

%% Distance between Tx and BD1 

for i = 1:length(Z)
    d1(i) = distance(t_x(1),Z(i,1),t_x(2),Z(i,2));
%     d1(i) = distance(t_x(1),Z(2,1),t_x(2),Z(2,2));
end

%% Distance between BD1 and Rx

for i = 1:length(Z)
    d2(i) = distance(r_x(1),Z(i,1),r_x(2),Z(i,2));
end

%% plotting concentric circles for different SINR values

% radii = SINR(200); %gathering the SINR values and utilising them as the ...
...radius for the circle
% angle = 0:0.001:2*pi; %utilising to draw the circle (different from theta)
% figure
% hold on
% for i = 1:length(radii)
%     xp = radii(i)*cos(angle);
%     yp = radii(i)*sin(angle);
%     plot(x+xp,y+yp,'Linewidth',1.5);
% end

%% 

% D1 = [0,x,30*lambda];
% D2 = [0,y,0];
% plot(D1,D2,'x-')

figure
scatter(t_x(1),t_x(2),'d');
hold on
scatter(r_x(1),r_x(2),'o');
plot(Z(:,1),Z(:,2),'rx')
b1 = [50,50];
b2 = [50,-50];
plot(b1,b2,'-')

xlabel('link distance (m)')
axis([-30 150 -50 50])
grid on

% scatter(5,30,'o');
% scatter(15,20,'o');
% scatter(15.77,23,'o');
% scatter(20,19,'o');
% scatter(15,36,'o');
% scatter(27,30,'o');
% scatter(15,15,'o');
% scatter(23,25,'o');


%% function for calculating the distance between any two points

function [d] = distance(x1,x2,y1,y2)
    d = sqrt((x2-x1)^2 + (y2-y1)^2);
end