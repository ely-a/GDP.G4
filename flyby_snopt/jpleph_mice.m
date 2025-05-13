function [r, v] = jpleph_mice(et, ntarg, ncent, frame)

% reads the jpl planetary ephemeris and gives the position and velocity
% of the point 'ntarg' with respect to point 'ncent' using MICE routines

% input

%   et    = TDB julian day at which interpolation is wanted

%   ntarg = integer number of 'target' celestial body

%   ncent = integer number of 'observer' celestial body

%   the numbering convention for 'ntarg' and 'ncent' is:

%        1 = mercury (199)      8 = neptune BC (8)
%        2 = venus (299)        9 = pluto BC (9)
%        3 = earth (399)       10 = moon (301)
%        4 = mars BC (4)       11 = sun (10)
%        5 = jupiter BC (5)    12 = solar system barycenter (0)
%        6 = saturn BC (6)
%        7 = uranus BC (7)     BC ==> barycenter

% output

%   r = position vector (kilometers)
%   v = velocity vector (kilometers/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set name of target body

switch ntarg

    case (1)

        targ = 'mercury';

    case (2)

        targ = 'venus';

    case (3)

        targ = 'earth';

    case (4)

        % Mars barycenter

        targ = '4';

    case (5)

        % Jupiter barycenter

        targ = '5';

    case (6)

        % Saturn barycenter

        targ = '6';

    case (7)

        % Uranus barycenter

        targ = '7';

    case (8)

        % Neptune barycenter

        targ = '8';

    case (9)

        % Pluto barycenter

        targ = '9';

    case (10)

        targ = 'moon';

    case (11)

        targ = 'sun';

    case (12)

        % solar system barycenter

        targ = '0';

end

% set name of central body

switch ncent

    case (1)

        obs = 'mercury';

    case (2)

        obs = 'venus';

    case (3)

        obs = 'earth';

    case (4)

        % Mars barycenter

        obs = '4';

    case (5)

        % Jupiter barycenter

        obs = '5';

    case (6)

        % Saturn barycenter

        obs = '6';

    case (7)

        % Uranus barycenter

        obs = '7';

    case (8)

        % Neptune barycenter

        obs = '8';

    case (9)

        % Pluto barycenter

        obs = '9';

    case (10)

        obs = 'moon';

    case (11)

        obs = 'sun';

    case (12)

        % solar system barycenter

        obs = '0';

end

% compute time, expressed as TDB seconds past J2000 TDB (2451545.0)

etime = 86400.0 * (et - 2451545.0);

% compute position and velocity vectors (no corrections)

starg = mice_spkezr(targ, etime, frame, 'NONE', obs);

% position vector (kilometers)

r = starg.state(1:3);

% velocity vector (kilometers/second)

v = starg.state(4:6);

