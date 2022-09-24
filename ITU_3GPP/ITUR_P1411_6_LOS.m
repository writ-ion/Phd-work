clc;
clear all;

%% Link A of backscatter communications

% Definition of parameters NLOS1

floors = 4;         % Floors in the building
h_r = 3*floors;   % height of the building (m)
orientation = 90;    % street orientation with respect to direct paths (degrees)
h_b = 30;            %height of base station (m)
h_m = 1.5;            %height of mobile station (m)
d = 70;              %distance from BS to MS (m)

% General parameters

f = [700 900 1800 2100 3500]*10^6;          %frequencies in Megahertz (MHz)
c = 3*10^8;                                 %speed of light
lambda = c./f;                               %wavelength (m)

% Path loss models

% LOS situation within street canyon (UHF)

L_bp = abs(20*log10(lambda.^2/(8*pi*h_b*h_m)));    %basic transmission loss 
R_bp = (4*h_b*h_m)./lambda;                        %breakpoint distance

L_los_l = ones(1,length(f));                        % Lower bound LOS
for i=1:length(f)
    if (d < R_bp(i) || d == R_bp(i))
        L_los_l(i) = L_bp(i) + 20*log10(d/R_bp(i));
    else
        L_los_l(i) = L_bp(i) + 40*log10(d/R_bp(i));
    end
end

L_los_u = ones(1,length(f));                        % Upper bound LOS
for i=1:length(f)
    if d > R_bp(i) 
        %|| d == R_bp(i)
        L_los_u(i) = L_bp(i) + 20 + 40*log10(d/R_bp(i));
    else
        L_los_u(i) = L_bp(i) + 20 + 25*log10(d/R_bp(i));
    end
end

L_los_m = ones(1,length(f));                        % middle bound LOS
for i=1:length(f)
    if d < R_bp(i) || d == R_bp(i)
        L_los_m(i) = L_bp(i) + 6 + 20*log10(d/R_bp(i));
    else
        L_los_m(i) = L_bp(i) + 6 + 40*log10(d/R_bp(i));
    end
end
% 
L_los_l
L_los_m
L_los_u