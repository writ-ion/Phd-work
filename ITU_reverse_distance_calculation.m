clc;
clear all;
close all;

L1 = 201.76;
L2 = 171.76;
f = 100;

E1 = 139.3-L1+20*log10(f);
E2 = 139.3-L2+20*log10(f);

x1 = (106.9-E1)/20;
x2 = (106.9-E2)/20;

d1 = 10^(x1)
d2 = 10^(x2)