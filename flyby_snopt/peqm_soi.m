function ydot = peqm_soi (t, y)

% flyby-planet-centered equations of motion

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global ip2 smu pmu jdate_ca

% current TDB julian date

jdate = jdate_ca + t / 86400.0;

% get current flyby planet heliocentric position vector

[rs2p, vs2p] = p2000(ip2, jdate);

% calculate position vector of sun relative to flyby planet

rp2s(1) = -rs2p(1);
rp2s(2) = -rs2p(2);
rp2s(3) = -rs2p(3);

rp2sm = norm(rp2s);

% compute heliocentric position vector of s/c

for i = 1:1:3
    
    rs2sc(i) = y(i) + rs2p(i);
    
end

rs2scm = norm(rs2sc);

rmag = norm(y(1:3));

r3 = rmag * rmag * rmag;

% acceleration due to the flyby planet

accp(1) = -pmu(ip2) * y(1) / r3; 

accp(2) = -pmu(ip2) * y(2) / r3;

accp(3) = -pmu(ip2) * y(3) / r3;

% acceleration due to the sun

accs(1) = -smu * (rs2sc(1) / rs2scm^3 + rp2s(1) / rp2sm^3);

accs(2) = -smu * (rs2sc(2) / rs2scm^3 + rp2s(2) / rp2sm^3);

accs(3) = -smu * (rs2sc(3) / rs2scm^3 + rp2s(3) / rp2sm^3);

% compute total flyby planet-centered acceleration vector

ydot(1) = y(4);

ydot(2) = y(5);

ydot(3) = y(6);

ydot(4) = accs(1) + accp(1);

ydot(5) = accs(2) + accp(2);

ydot(6) = accs(3) + accp(3);

ydot = ydot';



