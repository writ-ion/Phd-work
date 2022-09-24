function [Pi_self_total] = SINR_self_BD2()
f = 200;
lambda = (3*10^8)/(f*1e6);
SINR_lim = -10; % in dB

%% location of the BD1, Tx and Rx

t_x = [0 0];
r_x = [100*lambda 0];

% Different locations for the BD1

theta = 135;
d = (10*lambda);

BD1 = [(r_x(1) + d*cos(theta)) (r_x(2) + d*sin(theta))];% 85.5 0;116 6.5;115 1.5];
s = size (BD1);
s = s(1);

% Location of the BD2

% grid for BD2 around BD1

% x-axis limits
BD2_X1 = BD1(1)-(2*lambda); 
BD2_X2 = BD1(1)+(2*lambda);
% y-axis limits
BD2_Y1 = BD1(2)-(2*lambda);
BD2_Y2 = BD1(2)+(2*lambda);

% The grid about the x-y axis
BD2_X = BD2_X1:0.1:BD2_X2;
BD2_Y = BD2_Y1:0.1:BD2_Y2;
[BD2_XX,BD2_YY] = ndgrid(BD2_X,BD2_Y);
BD2Z = [BD2_XX(:),BD2_YY(:)];

% grid for BD2 around Rx

% x-axis limits
R_X1 = r_x(1)-(2*lambda); 
R_X2 = r_x(1)+(2*lambda);
% y-axis limits
R_Y1 = r_x(2)-(2*lambda);
R_Y2 = r_x(2)+(2*lambda);

% The grid about the x-y axis
R2_X = R_X1:0.1:R_X2;
R2_Y = R_Y1:0.1:R_Y2;
[R2_XX,R2_YY] = ndgrid(R2_X,R2_Y);
R2Z = [R2_XX(:),R2_YY(:)];

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

%% Distance between Tx and BD2 (with respect to BD1)

d3 = ones(1,length(BD2Z));
for i = 1:length(BD2Z)
    d3(i) = distance(t_x(1),BD2Z(i,1),t_x(2),BD2Z(i,2));
end

%% Distance between BD2 and Rx (with respect to BD1)

% d4 = ones(1,length(BD2Z));
for i = 1:length(BD2Z)
    d4(i) = distance(r_x(1),BD2Z(i,1),r_x(2),BD2Z(i,2));
end


%% Distance between BD1 and BD2 (with respect to BD1)

% bd2bd = ones(length(BD2Z),length(BD1));
j=1;
while j <= length(s)
    for i = 1:length(BD2Z)
        bd2bd(j,i) = distance(BD1(j,1),BD2Z(i,1),BD1(j,2),BD2Z(i,2));
    end
    j = j+1;
end

%% SINR and interference power calculation for all the BD2 points (wrt BD1)

SINR_CI = ones(1,length(BD2Z));
SINR_SI = ones(1,length(BD2Z));
SINR_TI = ones(1,length(BD2Z));
Pi_cross = ones(1,length(BD2Z));
Pi_self = ones(1,length(BD2Z));

j=1;
while j <= length(s)
for i = 1:length(BD2Z)
    [SINR_CI(j,i),SINR_SI(j,i),SINR_TI(j,i),Pi_cross(j,i),Pi_self(j,i)] = SINR_bd(f,d1(j),d2(j),d3(i),d4(i),bd2bd(j,i));
end
    j = j+1;
end

SINR_cross = reshape(SINR_CI(1,:),length(BD2_X),length(BD2_Y));
SINR_cross(SINR_cross < SINR_lim) = NaN;

SINR_self = reshape(SINR_SI(1,:),length(BD2_X),length(BD2_Y));
SINR_self(SINR_self < SINR_lim) = NaN;

Pi_self = sort(Pi_self);
Pi_self(end) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Distance between Tx and BD2 (with respect to Rx)

d3_r = ones(1,length(R2Z));
for i = 1:length(R2Z)
    d3_r(i) = distance(t_x(1),R2Z(i,1),t_x(2),R2Z(i,2));
end

%% Distance between BD2 and Rx (with respect to Rx)

% d4 = ones(1,length(BD2Z));
for i = 1:length(R2Z)
    d4_r(i) = distance(r_x(1),R2Z(i,1),r_x(2),R2Z(i,2));
end


%% Distance between BD1 and BD2 (with respect to Rx)

% bd2bd = ones(length(BD2Z),length(BD1));
j=1;
while j <= length(s)
    for i = 1:length(R2Z)
        bd2bd_r(j,i) = distance(BD1(j,1),R2Z(i,1),BD1(j,2),R2Z(i,2));
    end
    j = j+1;
end

%% SINR and interference power calculation for all the BD2 points (wrt Rx)

SINR_CI_r = ones(1,length(R2Z));
SINR_SI_r = ones(1,length(R2Z));
SINR_TI_r = ones(1,length(R2Z));
Pi_cross_r = ones(1,length(R2Z));
Pi_self_r = ones(1,length(R2Z));

j=1;
while j <= length(s)
for i = 1:length(R2Z)
    [SINR_CI_r(j,i),SINR_SI_r(j,i),SINR_TI_r(j,i),Pi_cross_r(j,i),Pi_self_r(j,i)] = SINR_bd(f,d1(j),d2(j),d3_r(i),d4_r(i),bd2bd_r(j,i));
end
    j = j+1;
end

SINR_self_r = reshape(SINR_SI_r(1,:),length(R2_X),length(R2_Y));
SINR_self_r(SINR_self_r < SINR_lim) = NaN;


Pi_self_total = [Pi_self,Pi_self_r];
% figure
% cdfplot(Pi_self)
% figure
% cdfplot(Pi_self_r)
% figure
% cdfplot(Pi_self_total)

%% Illustration of the environment
% CLIM = [0 60];
figure
imagesc(SINR_cross)
% mesh(BD2_XX,BD2_YY,SINR_cross);
% colorbar
a=colorbar;
ylabel(a,'SINR (dB)','FontSize',16,'Rotation',270);
xlabel('X [dm]')
ylabel('Y [dm]')
title('Cross Interference (at the BD1)')
% axis([0 60 0 60])

figure
imagesc(SINR_self)
% mesh(BD2_XX,BD2_YY,SINR_self);
% colorbar
a=colorbar;
ylabel(a,'SINR (dB)','FontSize',16,'Rotation',270);
xlabel('X [dm]')
ylabel('Y [dm]')
title('Self Interference (at the BD1)')

figure
imagesc(SINR_self_r)
% mesh(R2_XX,R2_YY,SINR_self_r);
% colorbar
a=colorbar;
ylabel(a,'SINR (dB)','FontSize',16,'Rotation',270);
xlabel('X [dm]')
ylabel('Y [dm]')
title('Self Interference (at the Rx)')
% axis([0 160 -50 50])
% grid on

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
% SINR_self_BD2R = [SINR_self;SINR_self_r];

% SI_cdf = reshape(Pi_self,[],1); %self interference
% CI_cdf = reshape(Pi_cross,[],1); %cross interference

% figure
% h(1, 1) = cdfplot(CI_cdf);
% hold on
% h(1, 2) = cdfplot(SI_cdf);
% set( h(:,1), 'LineStyle', ':', 'Color', 'r');
% set( h(:,2), 'LineStyle', '--', 'Color', 'b');
% xlabel('Interference Power (dBm)')
% legend ('Cross Interference','Self Interference')
% axis([0 0.5*10^-13 0.80 1])

% figure
% cdfplot(CI_cdf)
% xlabel('Cross Interference Power (dBm)')

% figure
% cdfplot(SI_cdf)
% xlabel('Self Interference Power (dBm)')

% toc
%% function for calculating the distance between any two points

function [d] = distance(x1,x2,y1,y2)
    d = sqrt((x2-x1)^2 + (y2-y1)^2);
end
end