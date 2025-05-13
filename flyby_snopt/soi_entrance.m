function [value, isterminal, direction] = soi_entrance(t, y)

% entrance to SOI of flyby planet event function

% input

%  t = time since departure (seconds)
%  y = spacecraft heliocentric state vector at time t 
%      (kilometers and kilometers/second)

% output

%  value = delta-SOI (kilometers)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ip2 rsoi jdtdb1 rp2sc vp2sc

% current TDB julian date

jdtdb = jdtdb1 + t / 86400.0;

% compute the heliocentric position and velocity of the flyby planet 
% (kilometers and kilometers/second)

[rfbp, vfbp] = p2000_ecl(ip2, jdtdb);

rsc = y(1:3);

vsc = y(4:6);

% state vector from flyby planet to spacecraft
% (kilometers and kilometers/second)

rp2sc = rfbp - rsc;

vp2sc = vfbp - vsc;

% objective function (delta-SOI distance; kilometers)

value = norm(rp2sc) - rsoi(ip2);

isterminal = 1;

direction =  [];


