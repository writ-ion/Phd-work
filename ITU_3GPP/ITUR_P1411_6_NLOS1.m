clc;
clear all;

%% Link A of backscatter communications

% Definition of parameters NLOS1

floors = 4;         % Floors in the building
h_r = 20;   % height of the building (m)
b = 20;             % average building separation (m)
w = 20;            % street width (m)
orientation = 90;    % street orientation with respect to direct paths (degrees)
h_b = 30;            %height of base station (m)
h_m = 1.5;            %height of mobile station (m)
l = 1000;              %length of path covered by buildings (m)
d = 1000;              %distance from BS to MS (m)

% General parameters

f = [700 900 1800 2100 3500]*10^6;          %frequencies in Megahertz (MHz)
c = 3*10^8;                                 %speed of light
lambda = c./f;                               %wavelength (m)


% NLOS Propagation over rooftops for urban area

f = [700 900 1800 2100 3500];

%FSPL

L_bf = 32.4 + 20*log10(d/1000) + 20*log10(f);       %free-space loss

%Rooftop to street

if orientation > -1 && orientation < 35
    l_ori = -10 + 0.354*orientation;
elseif orientation > 34 && orientation < 55
    l_ori = 2.5 + 0.075*(orientation - 35);
elseif orientation > 54 && orientation < 91
    l_ori = 4 - 0.114*(orientation - 55);
end

dh_m = h_r - h_m;

l_rts = -8.2 - 10*log10(w) + 10*log10(f) + 20*log10(dh_m) + l_ori;

%--------------------------------------------------------------------------------------------------------%
%Multiple screen diffraction

dh_b = h_b - h_r;
d_s = (lambda*(d^2))/dh_b^2;
d_bp = ones(1,length(f));

for i=1:length(f)
    d_bp(i) = abs(dh_b)*(sqrt(l/lambda(i))); %d_bp, lambda
end
% d_bp = floor(d_bp);

l1_msd = L1_MSD(d,dh_b,b,f,h_b,h_r);

l1_msd_bp = ones(1, length(f));
for i = 1:length(f)
    l1_msd_bp(i) = L1_MSD(d_bp(i),dh_b,b,f(i),h_b,h_r); %l1_msd_bp, d_bp, f
end

l2_msd = L2_MSD(d,dh_b,b,f,lambda,h_b,h_r);

l2_msd_bp = ones(1, length(f));
for i = 1:length(f)
    l2_msd_bp(i) = L2_MSD_BP(d_bp(i),dh_b,b,f(i),h_b,h_r); %l2_msd_bp, d_bp, f
end

%--------------------------------------------------------------------------------------------------------%

chi = 0.1;
upsilon = 0.0417;

for i = 1:length(f)
    l_upp(i) = l1_msd_bp(i);
    l_low(i) = l2_msd_bp(i);
    l_mid(i) = (l_upp(i) + l_low(i))/2;
    zeta(i) = (l_upp(i) - l_low(i))*upsilon;
    dh_bp(i) = l_upp(i) - l_low(i);
end

% Overall MSD

l_msd = length(f);
for i = 1:length(f)
    if l > d_s(i) && dh_bp(i) > 0
        l_msd(i) = -tanh((log(d)-log(d_bp(i)))/chi)*((l1_msd(i))-l_mid(i)) + l_mid(i);
    elseif (l < d_s(i) || l == d_s(i)) && dh_bp(i) > 0
        l_msd(i) = tanh((log(d)-log(d_bp(i)))/chi)*(l2_msd(i)-l_mid(i)) + l_mid(i);
    elseif dh_bp(i) == 0
        l_msd(i) = l2_msd(i);
    elseif l > d_s(i) && dh_bp(i) < 0
        l_msd(i) = l1_msd(i) - tanh((log(d)-log(d_bp(i)))/zeta(i))*((l_upp(i)-l_mid(i))) - l_upp(i) + l_mid(i);
    elseif (l < d_s(i) || l == d_s(i)) && dh_bp(i) < 0
        l_msd(i) = l2_msd(i) + tanh((log(d)-log(d_bp(i)))/zeta(i))*((l_mid(i)-l_low(i))) + l_mid(i) - l_low(i);
    end
end

L_NLOS1 = ones(1,length(f));
for i = 1:length(f)
    if l_rts + l_msd > 0
        L_NLOS1(i) = L_bf(i) + l_rts(i) + l_msd(i);
    else
        L_NLOS1(i) = L_bf(i);
    end
end
L_NLOS1
%--------------------------------------------------------------------------------------------------------%

%L1(msd) and L2(msd) definition:

%L1(msd)
function [l1_msd] = L1_MSD(d,dh_b,b,f,h_b,h_r)

%Lbsh definition
    if h_b > h_r
        l_bsh = -18*log10(1 + dh_b);
    else
        l_bsh = 0;
    end

%ka definition
    k_a = ones(1,length(f));
    for i=1:length(f)
        if h_b > h_r && f(i) > 2000 %f
            k_a(i) = 71.4; %k_a
        elseif (h_b < h_r || h_b == h_r) && f(i) > 2000 && (d > 500 || d == 500)
            k_a(i) = 73 - 0.8*dh_b;
        elseif (h_b < h_r || h_b == h_r) && f(i) > 2000 && d < 500
            k_a(i) = 73 - 1.6*dh_b*(d/1000);
        elseif h_b > h_r && (f(i) < 2000 || f(i) == 2000)
            k_a(i) = 54;
        elseif (h_b < h_r || h_b == h_r) && (f(i) < 2000 || f(i) == 2000) && (d > 500 || d == 500)
            k_a(i) = 54 - 0.8*dh_b;
        elseif (h_b < h_r || h_b == h_r) && (f(i) < 2000 || f(i) == 2000) && d < 500
            k_a(i) = 54 - 1.6*dh_b*(d/1000);
        end
    end

%kd definition
    if h_b > h_r
        k_d = 18;
    else
        k_d = 18 - 15*(dh_b/h_r);
    end

%kf definition
    k_f = ones(1,length(f));
    for i=1:length(f)
        if f(i) > 2000 %f
            k_f(i) = -8; %k_f
        else 
            k_f(i) = -4 + 0.7*((f(i)/925)-1); %k_f, f
        end
    end

    l1_msd = ones(1,length(f));
    for i=1:length(f)
        l1_msd(i) = l_bsh + k_a(i) + k_d*log10(d/1000) + k_f(i)*log10(f(i)) - 9*log10(b); %k_f, f
    end
end

%--------------------------------------------------------------------------------------------------------%

%L2(msd) 
function [l2_msd] = L2_MSD(d,dh_b,b,f,lambda,h_b,h_r)

    theta = atan(dh_b/b);
    rho = sqrt(dh_b^2 + b^2);

    dh_u = ones(1,length(f));
    for i=1:length(f)
        dh_u(i) = 10^(-log10(sqrt(b/lambda(i)))-(log10(d)/9)+((10/9)*log10(b/2.35))); %dh_u, lambda
    end

    dh_l = ones(1,length(f));
    for i=1:length(f)
        dh_l(i) = ((0.00023*b^2 - 0.1827*b - 9.4978)/(log10(f(i))^2.938)) + 0.000781*b + 0.06923; %dh_l, f
    end

    q_m = ones(1,length(f));
    for i=1:length(f)
        if h_b > (h_r + dh_u(i)) %dh_u
            q_m(i) = 2.35*((dh_b/d)*(sqrt(b/lambda(i))))^0.9; %q_m, lambda
        elseif (h_b < h_r + dh_u(i) || h_b == h_r + dh_u(i)) && (h_b > h_r + dh_l(i) || h_b == h_r + dh_l(i)) %dh_u, dh_l
            q_m(i) = b/d; %q_m
        elseif h_b < (h_r + dh_l(i)) %dh_l
            q_m(i) = (b/(2*pi*d)*(sqrt(lambda(i)/rho))*((1/theta)-(1/(2*pi + theta)))); %q_m, lambda
        end
    end

    l2_msd = ones(1,length(f));
    for i=1:length(f)
        l2_msd(i) = -10*log10(q_m(i)^2); %l2_msd, q_m
    end
end
%-------------------------------------------------------------------------------------------%
%L2(msd) 
function [l2_msd_bp] = L2_MSD_BP(d,dh_b,b,f,h_b,h_r)

    theta = atan(dh_b/b);
    rho = sqrt(dh_b^2 + b^2);
    lambda = (3*10^8)/(f*10^6);

%     dh_u = ones(1,length(f));
%     for i=1:length(f)
        dh_u = 10^(-log10(sqrt(b/lambda))-(log10(d)/9)+((10/9)*log10(b/2.35)));
%     end

%     dh_l = ones(1,length(f));
%     for i=1:length(f)
        dh_l = ((0.00023*b^2 - 0.1827*b - 9.4978)/(log10(f)^2.938)) + 0.000781*b + 0.06923;
%     end

%     q_m = ones(1,length(f));
%     for i=1:length(f)
        if h_b > (h_r + dh_u)
            q_m = 2.35*((dh_b/d)*(sqrt(b/lambda)))^0.9;
        elseif (h_b < h_r + dh_u || h_b == h_r + dh_u) && (h_b > h_r + dh_l || h_b == h_r + dh_l)
            q_m = b/d;
        elseif h_b < (h_r + dh_l)
            q_m = (b/(2*pi*d)*(sqrt(lambda/rho))*((1/theta)-(1/(2*pi + theta))));
        end
%     end

%     l2_msd_bp = ones(1,length(f));
%     for i=1:length(f)
        l2_msd_bp = -10*log10(q_m^2);
%     end
end