function [value, isterminal, direction] = fpa_event(~, y)

% flight path angle wrt flyby planet event function

% input

%  y = spacecraft planet-centered state vector 
%      (kilometer and kilometer/second)

% output

%  value = sine of flight path angle wrt flyby planet

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract the spacecraft position and velocity vectors
% relative to the flyby planet (kilometer and kilometer/second)

rsc = y(1:3);

vsc = y(4:6);

% sine of the flight path angle (radians)

value = dot(rsc, vsc) / (norm(rsc) * norm(vsc));

isterminal = 1;

direction =  [];


