% Device to device link
clc;
clear all;
close all;

d = 1:1:100;
% d = 1:1:99;
% f = [700 900 1800 2100 3500];
f = 700;
sigma = 7;
w = 20;
p = 1:1:99;
% p = 50;

% q = 1:1:99;
% p_los = [-11.3 -7.9 0.0 10.6 20.3];
% p_los_lin = 10.^(p_los/10);
% p_nlos = [-16.3 -9.0 0.0 9.0 16.3];
% p_nlos_lin = 10.^(p_nlos/10);
% d_los = [976 276 44 16 10];

%LOS loss (Step 1-3)

L_los_total = LOS_loss(f,d,p,sigma);

%NLOS loss (Step 4-6)

L_nlos_total = NLOS_loss(f,d,p,sigma);

%Step 7
for i = 1:99
    if p(i) < 45
        d_los_p(i) = 212*(log10(p(i)./100)).^2 - 64*log10(p(i)./100);
    else
        d_los_p(i) = 79.2 - 70*(p(i)./100);
    end
end

%Step 8
L_los = LOS_loss(f,d_los_p,p,sigma);
L_nlos = NLOS_loss(f,d_los_p + w,p,sigma);

L = ones(99);
for i = 1:length(d)
    for j = 1:length(p)
    if d(i) < d_los_p(j)
        L(j,i) = L_los_total(j,i);
    else%if d(i) > d_los_p(j) + w
        L(j,i) = L_nlos_total(j,i);
%     else
%         L(j,i) = L_los(j,i) + (L_nlos(j,i) - L_los(j,i))*(d(i) - d_los_p(j))/w;
    end
    end
end

i = [1 10 50 90 99];
plot(d,L(i,:));
set(gca, 'YDir','reverse')
grid on;
hold on;
legend('1%','10%','50%','90%','99%')
xlabel('Distance (m)')
ylabel('Total loss (dB)')


% plot(d,L(10,:));
% plot(d,L(99,:))
% L'
%L_los_i = interp1(p,d_los,q)

%LOS calculation function
function [L_los_final] = LOS_loss(f,d,p,sigma)
    
    L_los_d = 32.45 + 20*log10(f) + 20*log10(d./1000);

    %LOS location correction
    L_los_p = 1.5624*sigma*(sqrt(-2*log(1-p./100)) - 1.1774);

    L_los_final = ones(1,length(d));
    %Total LOS loss
    for i = 1:length(d)
        for j = 1:length(p)
        L_los_final(j,i) = L_los_d(i) + L_los_p(j);
        end
    end
%     L_los_final = L_los;
end

%NLOS calculation function
function [L_nlos_final] = NLOS_loss(f,d,p,sigma)
    
    %median value for NLOS loss
    L_nlos_d = 9.5 + 45*log10(f) + 40*log10(d./1000) + 6.8; %urban considered here

    %NLOS location correction
    L_nlos_p = sigma*norminv(p/100);

    L_nlos_final = ones(1,length(d));
    %Total NLOS loss
    for i = 1:length(d)
        for j = 1:length(p)
            L_nlos_final(j,i) = L_nlos_d(i) + L_nlos_p(j);
        end
    end
%     L_nlos_final = L_nlos;
end