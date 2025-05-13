function oev = flyby_hyper(mu, vinfi, vinfo, rp)

% classical orbital elements of a planet-centered
% flyby hyperbola at periapsis passage

% input

%  mu    = flyby planet gravitational constant (kilometers^3/seconds^2)
%  vinfi = incoming v-infinity vector (kilometers/second)
%  vinfo = outgoing v-infinity vector (kilometers/second)
%  rp    = planet-centered periapsis distance (kilometers)

% output

%  oev(1) = semimajor axis (kilometers)
%  oev(2) = orbital eccentricity (non-dimensional)
%  oev(3) = orbital inclination (radians)
%  oev(4) = argument of perigee (radians)
%  oev(5) = right ascension of ascending node (radians)
%  oev(6) = true anomaly (radians)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pv = zeros(3, 1);

% magnitude of incoming v-infinity vector (kilometers/second)

vinfim = norm(vinfi);

% semimajor axis (kilometers)

sma = -mu / (vinfim * vinfim);

% orbital eccentricity

ecc = 1.0 - rp / sma;

% unit vector in direction of asymptote

shat = vinfi / vinfim;

% unit vector normal to orbit plane of flyby hyperbola

wv = cross(vinfi, vinfo);

wvm = norm(wv);

what = wv / wvm;

% orbital inclination (radians)

orbinc = acos(what(3));

bv = cross(shat, what);

bvm = norm(bv);

bhat = bv / bvm;

% unit vector along z-axis of fundamental coordinate system

ak(1) = 0.0;
ak(2) = 0.0;
ak(3) = 1.0;

% unit vector in direction of ascending node

nv = cross(ak, what);

nvm = norm(nv);

nhat = nv / nvm;
 
% right ascension of the ascending node (radians)

raan = atan3(nhat(2), nhat(1));

% sine and cosine of true anomaly at infinity

cta = -1.0 / ecc;

sta = sqrt(1.0 - 1.0 / (ecc * ecc));

% unit vector in direction of periapsis

for i = 1:1:3
    
    pv(i) = -cta * shat(i) + sta * bhat(i);
    
end

pvm = norm(pv);

phat = pv / pvm;

% sine and cosine of argument of periapsis

cargper = dot(nhat, phat);

sargper = sqrt(1.0 - cargper * cargper);

if (phat(3) < 0.0)
    
   sargper = -sargper;
   
end

% argument of periapsis (radians)

argper = atan3(sargper, cargper);

% load orbital elements array

oev(1) = sma;
oev(2) = ecc;
oev(3) = orbinc;
oev(4) = argper;
oev(5) = raan;
oev(6) = 0.0;
