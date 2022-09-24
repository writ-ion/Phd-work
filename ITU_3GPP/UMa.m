% clc;
% clear all;
% close all;
function [L_total] = UMa(f,d)

% f = [700 900 1800 2100 3500];
% d = 30:3:327;
h_ut = 1.5;

if h_ut < 13 || h_ut == 13
    C_hut = 0;
else
    C_hut = ((h_ut - 13)/10)^1.5;
end

for i = 1:length(d)
    for j = 1:length(f)
        L_los(i,j) = 28.0 + 20*log10(f(j)/1000) + 22*log10(d(i));
        L_nlos(i,j) = 32.4 + 20*log10(f(j)/1000) + 30*log10(d(i));
        
        %LOS probability
        if d(i) < 18 || d(i) == 18
            p_los(i) = 1;
        else
            p_los(i) = (18/d(i) + exp(-d(i)/63)*(1-(18/d(i))))*(1 + C_hut*(5/4)*((d(i)/100)^3)*exp(-d(i)/150));
        end
        
        %Total loss after the probability is considered
        L_total(i,j) = p_los(i)*L_los(i,j) + (1-p_los(i))*L_nlos(i,j);
    end
end

% figure
% plot(d,L_los)
% grid on
% xlabel('distance (km)')
% ylabel('Path loss (dB)')
% title('UMi LOS path loss')
% legend('700','900','1800','2100','3500')
% 
% figure
% plot(d,L_nlos)
% grid on
% xlabel('distance (km)')
% ylabel('Path loss (dB)')
% title('UMi NLOS path loss')
% legend('700','900','1800','2100','3500')

% figure
% plot(d,L_total)
% grid on
% xlabel('distance (km)')
% ylabel('Path loss (dB)')
% title('UMa total path loss')
% legend('700','900','1800','2100','3500')
end