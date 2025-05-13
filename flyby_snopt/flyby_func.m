function [y, g] = flyby_func(x)

% patched-conic gravity assist mission constraints and objective function

% mean ecliptic and equinox of J2000 system

% input

%  x(1) = departure julian day - jdtdb0
%  x(2) = flyby julian day - jdtdb0
%  x(3) = arrival julian day - jdtdb0

% output

%  y(1) = delta-v objective function (kilometers/second)
%  y(2) = flyby altitude (kilometers)
%  y(3) = v-infinity matching (kilometers/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global smu ip1 ip2 ip3 jdtdb0 pmu req 

global dv_departure dvm_departure dv_arrival dvm_arrival

global flyby_sma flyby_ecc flyby_rp flyby_alt

global c3departure c3arrival rito1 vito1 rito2 vito2

global otype dvh vinfm_in vinfm_out

global vinf_in vinf_out rinf_in

y = zeros(3, 1);

sv1 = zeros(6, 1);

sv2 = zeros(6, 1);

% tdb Julian dates of current iteration

jdtdb1 =  jdtdb0 + x(1);

jdtdb2 =  jdtdb0 + x(2);

jdtdb3 =  jdtdb0 + x(3);

% posigrade transfer

direct = 1;

% less than one rev transfer

revmax = 0;

% compute heliocentric state vector of departure planet at departure 
% (ecliptic; kilometers and kilometers/second)

[rplanet1, vplanet1] = p2000_ecl(ip1, jdtdb1);

% compute heliocentric state vector of flyby planet at flyby
% (ecliptic; kilometers and kilometers/second)

[rplanet2, vplanet2] = p2000_ecl(ip2, jdtdb2);

% compute heliocentric state vector of arrival planet at arrival
% (ecliptic; kilometers and kilometers/second)

[rplanet3, vplanet3] = p2000_ecl(ip3, jdtdb3);

%--------------------------------------
% solve Lambert's problem for first leg
%--------------------------------------

tof1 = 86400.0 * (jdtdb2 - jdtdb1);
   
sv1(1:3) = rplanet1(1:3);
    
sv1(4:6) = vplanet1(1:3);
    
sv2(1:3) = rplanet2(1:3);
    
sv2(4:6) = vplanet2(1:3);

[vito1, vfto1] = glambert(smu, sv1, sv2, direct * tof1, revmax);

% load initial and final heliocentric position vectors for this leg
% (ecliptic; kilometers)

rito1 = rplanet1;

rfto1 = rplanet2;

% compute departure heliocentric delta-v vector, magnitude and energy

dv_departure = vito1 - vplanet1;

dvm_departure = norm(dv_departure);

c3departure = dvm_departure * dvm_departure;

% calculate incoming v-infinity vector and magnitude
% (ecliptic; kilometers/second)

vinf_in = vfto1 - vplanet2;

vinfm_in = norm(vinf_in);

rinf_in = rfto1 - rplanet2;

%---------------------------------------
% solve Lambert's problem for second leg
%---------------------------------------

tof2 = 86400.0 * (jdtdb3 - jdtdb2);
    
sv1(1:3) = rplanet2(1:3);

sv1(4:6) = vplanet2(1:3);

sv2(1:3) = rplanet3(1:3);

sv2(4:6) = vplanet3(1:3);

[vito2, vfto2] = glambert(smu, sv1, sv2, direct * tof2, revmax);

% load initial and final heliocentric position vectors for this leg
% (ecliptic; kilometers)

rito2 = rplanet2(1:3);

% heliocentric speed change due to flyby (kilometers/second)

dv = vito2 - vfto1;

dvh = norm(dv);

% calculate outgoing v-infinity vector and magnitude
% (ecliptic; kilometers/second)

vinf_out = vito2 - vplanet2;

vinfm_out = norm(vinf_out);

% determine flyby angle (radians)

v1crossv2 = cross(vinf_in, vinf_out);

v1crossv2m = norm(v1crossv2);

phi1 = 0.5 * pi + 0.5 * asin(v1crossv2m / (vinfm_in * vinfm_out));

% calculate flyby planet-centered hyperbolic orbital elements

flyby_ecc = -1.0 / cos(phi1);

flyby_sma = - pmu(ip2) / (vinfm_in * vinfm_in);

flyby_rp = flyby_sma * (1.0 - flyby_ecc);

% compute flyby altitude (kilometers)

flyby_alt = flyby_rp - req(ip2);

% compute arrival heliocentric delta-v vector and magnitude
% (ecliptic; kilometers/second)

dv_arrival = vplanet3 - vfto2;

dvm_arrival = norm(dv_arrival);

c3arrival = dvm_arrival * dvm_arrival;

% load scalar objective function

switch otype
    
case 1
    
   % departure delta-v (kilometer/second)
   
   f = dvm_departure;
   
case 2
    
   % arrival delta-v (kilometers/second)
   
   f = dvm_arrival;
   
case 3
    
   % departure + arrival delta-v (kilometers/second)
   
   f = dvm_departure + dvm_arrival;
   
end

% load objective function

y(1) = f;

% load mission constraints 
% (normalized flyby altitude and v-infinity matching)

y(2) = flyby_alt / req(ip2);

y(3) = vinfm_in - vinfm_out;

% no derivatives

g = [];
