clc 
clear
close all

% Sensitivity Analysis

%R1 = 149.6e6;

syms mu_S theta1 theta2 R1 v_per v_A dV_A R2 h

e = (R2-R1)/(R1*cosd(theta1) - R2*cosd(theta2));
v_r = (mu_S/h)*e*sind(theta1);
v_per = sqrt((v_A^2) - (v_r^2));
h_expr = v_per*R1;
h_expr = solve(h_expr,h);
h_expr = h_expr(1,1);

R2 = ((h^2)/mu_S)/(1+e*cosd(theta1));
R2_expr = subs(R2, h, h_expr);
solR2 = simplify( solve( R2 == R2_expr, R2 ) );
% 
% % delta R2 = (dR2/dV_A) *deltaV_A
% R2_dV_A =  diff(R2,v_A)*dV_A;



