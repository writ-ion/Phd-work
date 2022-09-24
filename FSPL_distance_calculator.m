clc;
clear all;
close all;

FSPL = 44; % in dB (path_loss = 77.8 - (-123.9752))
f=100; %in MHz
d=10^((FSPL-32.44-20*log10(f))/20) %in km

%d = 2.9293e+06