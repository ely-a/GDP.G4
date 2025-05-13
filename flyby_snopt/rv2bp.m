function [bplane, rv, tv, ibperr] = rv2bp(cbmu, r, v)

% convert planetocentric cartesian state vector 
% to b-plane elements

% input

%  cbmu = gravitational constant
%  r    = position vector
%  v    = velocity vector

% output

%  bplane(1) = vhp mag
%  bplane(2) = declination of asymptote
%  bplane(3) = right ascension of asymptote
%  bplane(4) = radius magnitude
%  bplane(5) = periapsis radius
%  bplane(6) = b-plane angle
%  bplane(7) = b dot t
%  bplane(8) = b dot r 
%  bplane(9) = b magnitude
%  bplane(10) = b-plane vector x-component
%  bplane(11) = b-plane vector y-component
%  bplane(12) = b-plane vector z-component
%  rv = b-plane r-vector
%  tv = b-plane t vector
%  ibperr = error flag (1 ==> not hyperbola)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ibperr = 0;

% angular momentum vector

wv = cross(r, v);

c1 = norm(wv);

wv = uvector(wv);

rrd = dot(r, v);

% position magnitude

rm = norm(r);

% velocity magnitude

vx = norm(v);

% hyperbolic speed

vhe2 = vx * vx - 2.0 * cbmu / rm;

if (vhe2 < 0.0)
   % orbit is not hyperbolic; set error flag
   
   ibperr = 1;
   
   tv = zeros(3);
   
   rv = zeros(3);
   
   bplane = zeros(12);
   
   return
end

vhem = sqrt(vhe2);

bplane(1) = vhem;

rdot = rrd / rm;

p = c1 * c1 / cbmu;

% semimajor axis

a = rm / (2.0 - rm * vx^2 / cbmu);

% orbital eccentricity

e = sqrt(1.0 - p / a);

% cosine and sine of true anomaly

cta = (p - rm) / (e * rm);

sta = rdot * c1 / (e * cbmu);

b = sqrt(p * abs(a));

ab = sqrt(a * a + b * b);

z = rm / c1 * v - rdot / c1 * r;

pv = cta * r / rm - sta * z;

qv = sta * r / rm + cta * z;
    
sv = -a / ab * pv + b / ab * qv;
    
bv = b * b / ab * pv + a * b / ab * qv;

sv = uvector(sv);

% declination of asymptote

bplane(2) = asin(sv(3));

% right ascension of asymptote

bplane(3) = atan3(sv(2), sv(1));

bplane(4) = rm;

% perapsis radius of hyperbola

rp = (e - 1.0) * cbmu / vhe2;

bplane(5) = rp;

ab = sqrt(sv(1) * sv(1) + sv(2) * sv(2));

% t vector

tv(1) = sv(2) / ab;

tv(2) = -sv(1) / ab;

tv(3) = 0.0;

rv = cross(sv, tv);

ab = norm(rv);

rv = uvector(rv);

% b dot t

bdott = dot(bv, tv);

bplane(7) = bdott;

% b dot r

bdotr = dot(bv, rv);

bplane(8) = bdotr;

% b-plane angle (radians)

theta = atan3(bdotr, bdott);

bplane(6) = theta;

% b magnitude

bplane(9) = rp * sqrt(1.0 + 2.0 * cbmu / (rp * vhe2));

% b-plane vector

bplane(10:12) = bv;


