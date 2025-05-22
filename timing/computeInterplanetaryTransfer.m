function result = computeInterplanetaryTransfer(departurePlanet, arrivalPlanet, launchDate, arrivalDate, plotFlag)
% computeInterplanetaryTransfer
% Solves the interplanetary Lambert problem and returns transfer orbit info
%
% Inputs:
%   departurePlanet - string (e.g. 'Earth')
%   arrivalPlanet   - string (e.g. 'Neptune')
%   launchDate      - datetime
%   arrivalDate     - datetime
%   plotFlag        - true/false to control plotting
% 
% Output:
%   result - struct with delta-v, orbital elements, and arrival info

    global mu
    mu = 1.32712440018e11; % km^3/s^2, Sun
    aunit = 149597870.691; % km

    % Planet name mapping
    planetList = {'Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune','Pluto'};
    ip1 = find(strcmpi(departurePlanet, planetList));
    ip2 = find(strcmpi(arrivalPlanet, planetList));

    if isempty(ip1) || isempty(ip2)
        error('Invalid planet name. Use one of: %s', strjoin(planetList, ', '));
    end

    % Julian Dates
    jdate1 = juliandate(launchDate);
    jdate2 = juliandate(arrivalDate);

    tof = (jdate2 - jdate1) * 86400; % seconds
    taud = jdate2 - jdate1;          % days

    % Get planetary state vectors
    [ri, vi] = planet(ip1, jdate1);
    [rf, vf] = planet(ip2, jdate2);

    % Lambert solution
    sv1 = [ri(:); vi(:)];
    sv2 = [rf(:); vf(:)];
    [vito, vfto] = lambert_gooding(mu, sv1, sv2, tof, 0);

    % Orbital elements
    oevi = eci2orb1(mu, ri, vi);
    oevf = eci2orb1(mu, rf, vf);
    oevtoi = eci2orb1(mu, ri, vito);
    oevtof = eci2orb1(mu, rf, vfto);

    % Delta-v
    dvi = vito(:) - vi(:);
    dvf = vf(:) - vfto(:);

    % Output struct
    result = struct();
    result.departureDate = launchDate;
    result.arrivalDate = arrivalDate;
    result.tof_days = taud;
    result.v1 = vito;
    result.v2 = vfto;
    result.deltaV1 = dvi;
    result.deltaV2 = dvf;
    result.oe_departurePlanet = oevi;
    result.oe_transfer_initial = oevtoi;
    result.oe_transfer_final = oevtof;
    result.oe_arrivalPlanet = oevf;
    result.dv1_mag = norm(dvi);
    result.dv2_mag = norm(dvf);
    result.energy = norm(dvi)^2;

    % Optional Plotting
    if plotFlag 
        plotInterplanetaryTransfer(mu, aunit, ip1, ip2, jdate1, jdate2, ri, vito, departurePlanet, arrivalPlanet);
    end
end
