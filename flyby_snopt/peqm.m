function ydot = peqm (t, y)

% first-order flyby-planet-centered equations of motion

% NOTE: includes point-mass gravity of the sun

% input

%  t = time since SOI (seconds)
%  y = planet-centered state vector at time t (kilometers and kilometers/second)

% output

%  ydot = integration vector (kilometers/second/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rs2sc = zeros(3, 1);

global ip2 smu pmu jdtdb_soi_in

% current TDB julian date

jdtdb = jdtdb_soi_in + t / 86400.0;

% get current flyby planet heliocentric state vector (kilometer and kilometer/second)

[rs2p, ~] = p2000(ip2, jdtdb);

% calculate position vector of sun relative to flyby planet (kilometers)

rp2s(1) = -rs2p(1);
rp2s(2) = -rs2p(2);
rp2s(3) = -rs2p(3);

rp2sm = norm(rp2s);

% compute heliocentric position vector of spacecraft (kilometers)

for i = 1:1:3
    
    rs2sc(i) = y(i) + rs2p(i);
    
end

rs2scm = norm(rs2sc);

rmag = norm(y(1:3));

r3 = rmag * rmag * rmag;

% acceleration due to the flyby planet (kilometers/second/second)

accp(1) = -pmu(ip2) * y(1) / r3; 

accp(2) = -pmu(ip2) * y(2) / r3;

accp(3) = -pmu(ip2) * y(3) / r3;

% acceleration due to the sun (km/sec^2)

accs(1) = -smu * (rs2sc(1) / rs2scm^3 + rp2s(1) / rp2sm^3);

accs(2) = -smu * (rs2sc(2) / rs2scm^3 + rp2s(2) / rp2sm^3);

accs(3) = -smu * (rs2sc(3) / rs2scm^3 + rp2s(3) / rp2sm^3);

% compute total flyby planet-centered acceleration vector (kilometers/second/second)

ydot(1) = y(4);

ydot(2) = y(5);

ydot(3) = y(6);

ydot(4) = accs(1) + accp(1);

ydot(5) = accs(2) + accp(2);

ydot(6) = accs(3) + accp(3);

ydot = ydot';



