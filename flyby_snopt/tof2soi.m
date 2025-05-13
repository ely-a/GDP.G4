function dtof = tof2(mu, req, sma, ecc, r1, r2)

% time of flight between two radii

% input

%  sma = semimajor axis (kilometers)
%  ecc = eccentricity (non-dimensional)
%  r1  = initial radius (kilometers)
%  r2  = final radius (kilometers)

% output

%  dtof = time-of-flight between r1 and r2 (seconds)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pi2 = 2 * pi;

% orbital period

tau = pi2 * sqrt(sma * sma * sma / mu);

% time of flight

dtof1 = sqrt(sma / mu) * (sma * acos((sma - r1) / (sma * ecc)) ...
   - sqrt(2 * sma * r1 - sma * sma * (1 - ecc * ecc) - r1 * r1));

dtof2 = sqrt(sma / mu) * (sma * acos((sma - r2) / (sma * ecc)) ...
   - sqrt(2 * sma * r2 - sma * sma * (1 - ecc * ecc) - r2 * r2));
 
dtof = dtof2 - dtof1;

% if necessary, make positive

if (dtof < 0)
   dtof = dtof + tau;
end
