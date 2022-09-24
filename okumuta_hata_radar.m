clc;
clear all;
close all;

%% large city

f = input('Enter the frequency (in MHz): '); %fm frequency MHz
P_dbm = 46; %in dBm
G_t = 12; %gain of the antenna in dBi
A = 69.55;
B = 26.16;

h_t = input('Enter transmitter height (in meters): '); %transmitter height (meters)
h_r = input('Enter the receiver height: ');

gamma = input('Enter the path-loss exponent: ');
C = gamma*10 + 6.55*log10(h_t)
x = C - 6.55*log10(h_t)

Cm = input('Enter area correction factor: '); %area correction factor

condition = true;
% flag = 1;

while condition
    d = input('Enter distance between transmitter & receiver (in km): ');   %km distance from the transmitter/receiver to the sensor
    if isempty(d)
        condition = false;
        break;
    end
    
    a_hms1=3.2*(log10(11.75*h_r)).^2-4.97;
    
    L = A + B*log10(f) - 13.82*log10(h_t) - a_hms1 + (C - 6.55*log10(h_t))*log10(d) + Cm;
    %received power
    R1 = P_dbm + G_t - L;
%     figure
%     plot(R1,h_t,'-o');
%     legend('show')
%     ylabel("Height of transmitter (m)")
%     xlabel("Received signal strength (dBm)")
%     title("Height vs Received signal strength (large city)")
%     grid on;
%     hold on;
    
%     if flag == 1
%         figure
%         flag = 0;
%     end
    plot(L, h_t, '-x');
    legend('show')
    ylabel("Height of transmitter (m)")
    xlabel("Loss (dB)")
    title("Height vs Loss (large city)")
    hold on;
    grid on;
end
%% small city

% f = 100; %fm frequency MHz
% % Tx_pow = 60000; %60kW
% % P_dbm = 10*log10(Tx_pow)+30;
% P_dbm = 46; %in dBm
% G_t = 12; %gain of the antenna in dBi
% h_t = input('Enter transmitter/receiver height (in meters): '); %transmitter height (meters)
% h_r = input('Enter the transmitter/receiver height: ');
% 
% gamma = input('Enter the path-loss exponent: ');
% A = 69.55;
% B = 26.16;
% C = gamma*10 + 6.55*log10(h_t);
% Cm = input('Enter area correction factor: '); %-10 area correction factor
% condition = true;
% flag = 1;
% 
% while condition
%     d = input('Enter distance between transmitter & receiver (in km): '); %km distance from the transmitter/receiver to the sensor
%     if isempty(d)
%         condition = false;
%         break;
%     end
%     a_hms2=(1.1*log10(f)-0.7)*h_r - (1.56*log10(f)-0.8);
%     L2=A+B*log10(f)-13.82*log10(h_t)-a_hms2+(C-6.55*log10(h_t))*log10(d)+Cm;
%     %received power
%     R2=P_dbm+0-L2;%no gains considered
%     
%     plot(L2,h_t,'-x');
%     legend('show')
%     ylabel("Height of transmitter (m)")
%     xlabel("Loss (dB)")
%     title("Height vs Loss (small city)")
%     grid on;
%     hold on;
%     
% %     plot(R2,h_t,'-o');
% %     legend('show')
% %     ylabel("Height of transmitter (m)")
% %     xlabel("Received signal strength (dBm)")
% %     title("Height vs Received signal strength (small city)")
% %     grid on;
% %     hold on;
% end