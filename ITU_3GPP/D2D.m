% Device to device link
function [d2d_loss] = D2D(f,d)

% d = 1:1:100;
% d = 1:1:99;
% f = [700 900 1800 2100 3500];
% f = 700;
sigma = 7; 
w = 20;
% p = 1:1:99;
p = 50;

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

if p < 45
    d_los_p = 212*(log10(p/100))^2 - 64*log10(p/100);
else
    d_los_p = 79.2 - 70*(p/100);
end

%Step 8

L_los = LOS_loss(f,d_los_p,p,sigma);
L_nlos = NLOS_loss(f,d_los_p + w,p,sigma);

% L = ones(5,100);
for i = 1:length(d)
    for j = 1:length(f)
        if d(i) < d_los_p
            L(i,j) = L_los_total(j,i);
        else%if d(i) > d_los_p(j) + w
            L(i,j) = L_nlos_total(j,i);
%     else
%         L(j,i) = L_los(j,i) + (L_nlos(j,i) - L_los(j,i))*(d(i) - d_los_p(j))/w;
        end
    end
end

d2d_loss = L';

%figure
%plot(d,d2d_loss)
%grid on
%xlabel('distance (km)')
%ylabel('Path loss (dB)')
%title('ITU D2D path loss')
% legend('700','900','1800','2100','3500')


% i = [1 10 50 90 99];
% plot(d,L(i,:));
% set(gca, 'YDir','reverse')
% grid on;
% hold on;
% xlabel('Distance (m)')
% ylabel('Total loss (dB)')


% plot(d,L(10,:));
% plot(d,L(99,:))
% L'
%L_los_i = interp1(p,d_los,q)

%LOS calculation function
function [L_los_final] = LOS_loss(f,d,p,sigma)
    
    for i = 1:length(d)
        for j=1:length(f)
            L_los_d(j,i) = 32.45 + 20*log10(f(j)) + 20*log10(d(i)/1000);
        end
    end
    
    %LOS location correction
    L_los_p = 1.5624*sigma*(sqrt(-2*log(1-p./100)) - 1.1774);
    
    L_los_final = L_los_d + L_los_p;

%     L_los_final = ones(1,length(d));
%     %Total LOS loss
%     for i = 1:length(d)
%         for j = 1:length(f)
%             L_los_final(j,i) = L_los_d(i) + L_los_p;
%         end
%     end
%     L_los_final = L_los;
end

%NLOS calculation function
function [L_nlos_final] = NLOS_loss(f,d,p,sigma)
    
    %median value for NLOS loss
    for i=1:length(d)
        for j=1:length(f)
            L_nlos_d(j,i) = 9.5 + 45*log10(f(j)) + 40*log10(d(i)/1000) + 6.8; %urban considered here
        end
    end
    %NLOS location correction
    L_nlos_p = sigma*norminv(p/100);
    
    L_nlos_final = L_nlos_d + L_nlos_p;

%     L_nlos_final = ones(1,length(d));
%     %Total NLOS loss
%     for i = 1:length(d)
%         for j = 1:length(p)
%             L_nlos_final(j,i) = L_nlos_d(i) + L_nlos_p(j);
%         end
%     end
%     L_nlos_final = L_nlos;
end
end