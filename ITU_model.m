%% ITU model for Tx to sensor
clc;
clear all;
close all;
 
%% Simulation parameters
 
f=100; %MHz
d=30;  %distance to sensor
h1=248; %height of base station
h2=1; %height of sensor
 
%% Steps for application of the recommendation
 
%Step 1: Land path
%Step 2: 50% time
%Step 3: Frequency (100 MHz)
%Step 4: Distance (30 kms to sensor)
%Step 5: Continue
%Step 6: Continue
%Step 7: Continue
%Step 8:
    %Step 8.1: Antenna height greater than 10 meters
    %Step 8.1.1: Required h1 is 248 which does not coincide with any of the
    %nominal values
        h_inf=150; %effective height below h1
        h_sup=300; %effective height above h1
        E_inf=30; %dB values from graph
        E_sup=37; %dB
    %Step 8.1.2: continue
    %Step 8.1.3: continue
    %Step 8.1.4: continue
    %Step 8.1.5: required distance of 30 coincides with nominal valuesin
    %Table 1
    %Step 8.1.6:
        E=E_inf+(E_sup-E_inf)*(log10(h1/h_inf)/log10(h_sup/h_inf)); %Annex 5, 4.1
    %Step 8.2: antenna height above 10m
%Step 9: frequency of 100 MHz coincides
%Step 10: percentage time is 50%
%Step 11: only land path is considered
%Step 12: Terrain information not considered
%Step 13: Tropospheric scattering not considered (information unavailable)
%Step 14: Sensors located at height of 1m, have to be calculated
%Step 15: Clutter not taken into account as open aarea is considered
%Step 16: Slope path correction
    d_slope=sqrt(d^2+10^-6*(h1-h2)^2); %terrain information unavailable
    correction=10*log10(d/d_slope);
%Step 17: Path length greater than 1 km
%Step 18: Terrain information unavailable
%Step 19:
    %E_max=E_fs %for land paths, the max E value should not exceed this
    E_fs=106.9-20*log10(d);
%Step 20: Equivalent basic transmission loss
    L=139.3-E+20*log10(f);