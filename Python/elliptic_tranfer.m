clear 
close all
clc

R_earth = 6378;
mu_earth = 398600.4419;
r_H = 1.5e6;

alt_LEO = 300;
r_LEO = R_earth + alt_LEO;
v_LEO = sqrt(mu_earth/r_LEO);

r_p = r_LEO;
e = 0.9;
a = r_p/(1-e);
P_e = (2 * pi * sqrt(a^3/mu_earth))/(3600 * 24);
r_a = a*(1+e);
h = sqrt(mu_earth * a * (1 - e^2));
v_p = h/r_p;

dv_US = v_p - v_LEO;
dv_KS = 9.5 - dv_US;

v_a = h/r_a;
delta = 30; 
dv_plane = sqrt(2*v_a^2 - 2*v_a^2*cosd(delta));