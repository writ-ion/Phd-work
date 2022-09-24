% clc;
% clear all;
% close all;
function [L_NLOS1] = NLOS1(f,d)
% Link A of backscatter communications

% Definition of parameters NLOS1

% floors = 4;         % Floors in the building
h_r = 15;   % height of the building (m)
b = 20;             % average building separation (m)
w = 20;            % street width (m)
orientation = 30;    % street orientation with respect to direct paths (degrees)
h_b = 30;            %height of base station (m)
h_m = 1.5;            %height of mobile station (m)
l = 50;              %length of path covered by buildings (m)
% d = 30:3:330;              %distance from BS to MS (m)

% General parameters

% f = [700 900 1800 2100 3500]*10^6;          %frequencies in Megahertz (MHz)
% f = 700;
c = 3*10^8;                                 %speed of light
lambda = c./(f*10^6);                               %wavelength (m)


% NLOS Propagation over rooftops for urban area

% f = [700 900 1800 2100 3500];
% f = 700;

%FSPL
for i = 1:length(f)
    for j = 1:length(d)
        L_bf(j,i) = 32.4 + 20*log10(d(j)/1000) + 20*log10(f(i));       %free-space loss
    end
end

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

for i = 1:length(f)
    for j = 1:length(d)
        d_s(j,i) = (lambda(i)*(d(j)^2))/dh_b^2;
    end
end

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

l_msd = ones(length(d),length(f));
for i = 1:length(f)
    if l > d_s(j,i) && dh_bp(i) > 0
        l_msd(j,i) = -tanh((log(d(j))-log(d_bp(i)))/chi)*((l1_msd(j,i))-l_mid(i)) + l_mid(i);
    elseif (l < d_s(j,i) || l == d_s(j,i)) && dh_bp(i) > 0
        l_msd(j,i) = tanh((log(d(j))-log(d_bp(i)))/chi)*(l2_msd(j,i)-l_mid(i)) + l_mid(i);
    elseif dh_bp(i) == 0
        l_msd(j,i) = l2_msd(j,i);
    elseif l > d_s(j,i) && dh_bp(i) < 0
        l_msd(j,i) = l1_msd(j,i) - tanh((log(d(j))-log(d_bp(i)))/zeta(i))*((l_upp(i)-l_mid(i))) - l_upp(i) + l_mid(i);
    elseif (l < d_s(j,i) || l == d_s(j,i)) && dh_bp(i) < 0
        l_msd(j,i) = l2_msd(j,i) + tanh((log(d(j))-log(d_bp(i)))/zeta(i))*((l_mid(i)-l_low(i))) + l_mid(i) - l_low(i);
    end
end

L_NLOS1 = ones(length(d),length(f));
for i = 1:length(f)
    for j=1:length(d)
        if l_rts + l_msd > 0
            L_NLOS1(j,i) = L_bf(j,i) + l_rts(i) + l_msd(j,i);
        else
            L_NLOS1(j,i) = L_bf(j,i);
        end
    end
end

% plot(d,L_NLOS1)
% grid on;
% legend('700','900','1800','2100','3500')
% axis([0 99 60 130])
% title('Path loss ITU macro')
% xlabel('distance (m)')
% ylabel('Path loss (dB)')

% additional_loss = 0; %dB
% L_total = L_NLOS1+additional_loss;
% figure;
% plot(d,L_total)
% grid on;
% legend('700','900','1800','2100','3500')
% % axis([0 99 100 170])
% title('Path loss + additional loss')
% xlabel('distance (m)')
% ylabel('Path loss (dB)')
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
    k_a = ones(length(d),length(f));
    for i=1:length(f)
        for j=1:length(d)
            if h_b > h_r && f(i) > 2000 %f
                k_a(j,i) = 71.4; %k_a
            elseif (h_b < h_r || h_b == h_r) && f(i) > 2000 && (d(j) > 500 || d(j) == 500)
                k_a(j,i) = 73 - 0.8*dh_b;
            elseif (h_b < h_r || h_b == h_r) && f(i) > 2000 && d(j) < 500
                k_a(j,i) = 73 - 1.6*dh_b*(d/1000);
            elseif h_b > h_r && (f(i) < 2000 || f(i) == 2000)
                k_a(j,i) = 54;
            elseif (h_b < h_r || h_b == h_r) && (f(i) < 2000 || f(i) == 2000) && (d(j) > 500 || d(j) == 500)
                k_a(j,i) = 54 - 0.8*dh_b;
            elseif (h_b < h_r || h_b == h_r) && (f(i) < 2000 || f(i) == 2000) && d(j) < 500
                k_a(j,i) = 54 - 1.6*dh_b*(d(j)/1000);
            end
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

    l1_msd = ones(length(d),length(f));
    for i=1:length(f)
        for j=1:length(d)
            l1_msd(j,i) = l_bsh + k_a(j,i) + k_d*log10(d(j)/1000) + k_f(i)*log10(f(i)) - 9*log10(b); %k_f, f
        end
    end
end

%--------------------------------------------------------------------------------------------------------%

%L2(msd) 
function [l2_msd] = L2_MSD(d,dh_b,b,f,lambda,h_b,h_r)

    theta = atan(dh_b/b);
    rho = sqrt(dh_b^2 + b^2);

    dh_u = ones(length(d),length(f));
    for i=1:length(f)
        for j=1:length(d)
            dh_u(j,i) = 10^(-log10(sqrt(b/lambda(i)))-(log10(d(j))/9)+((10/9)*log10(b/2.35))); %dh_u, lambda
        end
    end

    dh_l = ones(1,length(f));
    for i=1:length(f)
        dh_l(i) = ((0.00023*b^2 - 0.1827*b - 9.4978)/(log10(f(i))^2.938)) + 0.000781*b + 0.06923; %dh_l, f
    end

    q_m = ones(1,length(f));
    for i=1:length(f)
        for j=1:length(d)
            if h_b > (h_r + dh_u(i)) %dh_u
                q_m(j,i) = 2.35*((dh_b/d(j))*(sqrt(b/lambda(i))))^0.9; %q_m, lambda
            elseif (h_b < h_r + dh_u(i) || h_b == h_r + dh_u(i)) && (h_b > h_r + dh_l(i) || h_b == h_r + dh_l(i)) %dh_u, dh_l
                q_m(j,i) = b/d(j); %q_m
            elseif h_b < (h_r + dh_l(i)) %dh_l
                q_m(j,i) = (b/(2*pi*d(j))*(sqrt(lambda(i)/rho))*((1/theta)-(1/(2*pi + theta)))); %q_m, lambda
            end
        end
    end

    l2_msd = ones(1,length(f));
    for i=1:length(f)
        for j=1:length(d)
            l2_msd(j,i) = -10*log10(q_m(j,i)^2); %l2_msd, q_m
        end
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
end