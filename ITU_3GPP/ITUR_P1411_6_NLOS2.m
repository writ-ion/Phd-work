clc;
clear all;

% Definition of parameters NLOS2

h_b = 30;            %height of base station (m)
h_m = 1.5;           %height of mobile station (m)
w1 = 20;             %street width at the position of the BS (m)
w2 = 20;             %street width at the position of the MS (m)
x1 = 30;             %distance BS to street crossing (m)
x2 = 30;             %distance MS to street crossing (m)
d = 1000;

% General parameters

f1 = [700 900 1800]*10^6;          %frequencies in Megahertz (MHz)
c = 3*10^8;                                 %speed of light
lambda1 = c./f1;                               %wavelength (m)

% Propagation within street canyons for 800-2000 MHz

% Reflection path loss

alpha = pi/3;
f_alpha = 3.86/(alpha^3.5);
l_r = 20*log10(x1 + x2) + x1*x2*(f_alpha/(w1*w2)) + 20*log10((4*pi)./lambda1);

% Diffraction path loss

d_a = (40/2*pi) * (atan(x2/w2) + atan(x1/w1) - pi/2);
l_d = 10*log10(x1*x2*(x1+x2)) + 2*d_a - 0.1*(90 - (alpha*(180/pi))) + 20*log10((4*pi)./lambda1);

L_NLOS2_1 = -10*log10((10.^(-l_d/10))+(10.^(-l_r/10)));

% Propagation within strret canyons for 2 to 16 GHz

f2 = [2100 3500]*10^6;
lambda2 = c./f2;
d_corner = 30;      %meters
l_corner = 20;      %dB
beta = 6;

if x2 > w1/2 +1 || (x2 < w1/2 + 1 + d_corner ||  x2 == w1/2 + 1 + d_corner)
    l_c = (l_corner/log10(1 + d_corner)) * log10(x2 - w1/2);
elseif x2 > w1/2 + 1 + d_corner
    l_c = l_corner;
end

if x2 > w1/2 + 1 + d_corner
    l_att = 10 * beta * log10((x1 + x2)/(x1 + w1/2 + d_corner));
else
    l_att = 0;
end

[lower,upper,middle] = LOS_SC(h_b,h_m,d,lambda2);

L_NLOS2_2_lower = lower + l_c + l_att;
L_NLOS2_2_upper = upper + l_c + l_att;
L_NLOS2_2_middle = middle + l_c + l_att;

[L_NLOS2_1 L_NLOS2_2_middle]

% LOS situation within street canyon (UHF)

function [L_los_l, L_los_u, L_los_m] = LOS_SC(h_b,h_m,d,lambda)
    L_bp = abs(20*log10(lambda.^2/(8*pi*h_b*h_m)));    %basic transmission loss at breakpoint
    R_bp = (4*h_b*h_m)./lambda;                        %breakpoint distance

    L_los_l = ones(1,length(lambda));                        % Lower bound LOS
    for i=1:length(lambda)
        if d < R_bp(i) || d == R_bp(i)
            L_los_l(i) = L_bp(i) + 20*log10(d/R_bp(i));
        else
            L_los_l(i) = L_bp(i) + 40*log10(d/R_bp(i));
        end
    end

    L_los_u = ones(1,length(lambda));                        % Upper bound LOS
    for i=1:length(lambda)
        if d < R_bp(i) || d == R_bp(i)
            L_los_u(i) = L_bp(i) + 20 + 25*log10(d/R_bp(i));
        else
            L_los_u(i) = L_bp(i) + 20 + 40*log10(d/R_bp(i));
        end
    end

    L_los_m = ones(1,length(lambda));                        % Lower bound LOS
    for i=1:length(lambda)
        if d < R_bp(i) || d == R_bp(i)
            L_los_m(i) = L_bp(i) + 6 + 20*log10(d/R_bp(i));
        else
            L_los_m(i) = L_bp(i) + 6 + 40*log10(d/R_bp(i));
        end
    end
end