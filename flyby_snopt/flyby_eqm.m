function ydot = flyby_eqm (t, y)

% first-order flyby-planet-centered equations of motion

% NOTE: includes point-mass gravity of the sun

% input

%  t = time since periapsis passage (seconds)
%  y = planet-centered state vector at periapsis
%      (kilometers and kilometers/second)

% output

%  ydot = integration vector (kilometers/second/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ip2 smu pmu jdtdb2

rs2sc = zeros(3, 1);

% current TDB julian day

jdtdb = jdtdb2 + t / 86400.0;

% flyby planet heliocentric state vector (kilometer and kilometer/second)

[rs2p, ~] = p2000_ecl(ip2, jdtdb);

% position vector of sun relative to flyby planet (kilometers)

rp2s = -rs2p;

rp2sm = norm(rp2s);

% compute heliocentric position vector of spacecraft (kilometers)

for i = 1:1:3

    rs2sc(i) = y(i) + rs2p(i);

end

rs2scm = norm(rs2sc);

rmag = norm(y(1:3));

r3 = rmag * rmag * rmag;

% acceleration due to the flyby planet (kilometers/second/second)

accp(1:3) = -pmu(ip2) * y(1:3) / r3;

% acceleration due to the sun (kilometers/second/second)

accs(1:3) = -smu * (rs2sc(1:3) / rs2scm^3 + rp2s(1:3) / rp2sm^3);

% total flyby planet-centered acceleration vector
% (kilometers/second/second)

ydot(1:3) = y(4:6);

ydot(4:6) = accs(1:3) + accp(1:3);

ydot = ydot';



