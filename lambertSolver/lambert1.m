% lambert1.m     January 26, 2024

% solution of the two-body Earth orbit Lambert problem

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

global mu

% Earth gravitational constant (kilometers^3/seconds^2)

mu = 398600.4415;

clc; home;

fprintf('\n           program lambert1\n');

fprintf('\n  < two-body Earth orbit lambert problem >\n');

% request orbital elements of initial orbit

fprintf('\nclassical orbital elements of the initial orbit\n');

oev1 = getoe([1;1;1;1;1;1]);

% request orbital elements of final orbit

fprintf('\nclassical orbital elements of the final orbit \n');

oev4 = getoe([1;1;1;1;1;1]);

% transfer time

while(1)

    fprintf('\nplease input the transfer time in minutes\n');

    ttmins = input('? ');

    if (ttmins > 0)

        break;

    end

end

% time of flight (seconds)

tof = 60.0 * ttmins;

% request transfer direction

while(1)

    fprintf('\n orbital direction\n');

    fprintf('\n  <1> posigrade\n');

    fprintf('\n  <2> retrograde');

    fprintf('\n\n selection (1 or 2)\n');

    direct = input('? ');

    if (direct == 1 || direct == 2)

        if (direct == 2)

            direct = -1;

        end

        break;

    end

end

% request number of orbits around the Earth

while(1)

    fprintf('\nplease input the maximum number of transfer orbits around the Earth\n');

    revmax = input('? ');

    if (revmax >= 0)

        break;

    end

end

% compute state vectors of initial and final orbits

[ri, vi] = orb2eci(mu, oev1);

[rf, vf] = orb2eci(mu, oev4);

% solve Lambert's problem using Gooding's method

sv1 = zeros(6, 1);

sv2 = zeros(6, 1);

sv1(1:3) = ri(1:3);

sv1(4:6) = vi(1:3);

sv2(1:3) = rf(1:3);

sv2(4:6) = vf(1:3);

[vito, vfto] = lambert_gooding(mu, sv1, sv2, direct * tof, revmax);

% orbital elements of the transfer orbit at each delta-v

oev2 = eci2orb1(mu, ri, vito);

oev3 = eci2orb1(mu, rf, vfto);

% delta-v vectors (kilometers/second)

dvi(1) = vito(1) - vi(1);
dvi(2) = vito(2) - vi(2);
dvi(3) = vito(3) - vi(3);

dvf(1) = vf(1) - vfto(1);
dvf(2) = vf(2) - vfto(2);
dvf(3) = vf(3) - vfto(3);

% print results

fprintf('\n         program lambert1\n');

fprintf('\n< two-body Earth orbit lambert problem >\n');

fprintf('\norbital elements of the initial orbit\n');

oeprint1(mu, oev1);

fprintf('\norbital elements of the transfer orbit after the initial delta-v\n');

oeprint1(mu, oev2);

fprintf('\norbital elements of the transfer orbit prior to the final delta-v\n');

oeprint1(mu, oev3);

fprintf('\norbital elements of the final orbit \n');

oeprint1(mu, oev4);

fprintf('\ninitial delta-v vector and magnitude\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * dvi(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * dvi(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * dvi(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second', 1000.0 * norm(dvi));

fprintf('\n\nfinal delta-v vector and magnitude\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * dvf(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * dvf(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * dvf(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second', 1000.0 * norm(dvf));

fprintf('\n\ntotal delta-v               %12.6f  meters/second', 1000.0 * (norm(dvi) + norm(dvf)));

fprintf('\n\ntransfer time               %12.6f  minutes \n\n', ttmins);

% perform primer vector initialization

[pvi, pvdi] = pv_iniz(mu, tof, ri, vito, dvi, dvf);

% number of graphic data points

npts = 300;

% plot behavior of primer vector magnitude

dt = tof / npts;

xp1 = zeros(npts + 1, 1);

yp1 = zeros(npts + 1, 1);

yp2 = zeros(npts + 1, 1);

for i = 1:1:npts + 1

    t = (i - 1) * dt;

    if (t == 0.0)

        % initial value of primer magnitude and derivative

        pvm = norm(pvi);

        pvdm = dot(pvi, pvdi) / pvm;

    else

        % primer vector and derivative magnitudes at time t

        [pvm, pvdm] = pvector(mu, ri, vito, pvi, pvdi, t);

    end

    % load data array

    xp1(i) = t / 60.0;

    yp1(i) = pvm;

    yp2(i) = pvdm;

end

figure(1);

hold on;

plot(xp1, yp1, '-r');

plot(xp1, yp1, '.r');

title('Primer Vector Analysis', 'FontSize', 16);

xlabel('simulation time (minutes)', 'FontSize', 12);

ylabel('primer vector magnitude', 'FontSize', 12);

grid;

% create primer graphics disk files

print('primer.eps', '-depsc');

print('primer.jpg', '-djpeg');

print('primer.tif', '-dtiff');

% plot behavior of magnitude of primer derivative

figure(2);

hold on;

plot(xp1, yp2, '-r');

plot(xp1, yp2, '.r');

title('Primer Vector Analysis', 'FontSize', 16);

xlabel('simulation time (minutes)', 'FontSize', 12);

ylabel('primer derivative magnitude', 'FontSize', 12);

grid;

% create primer derivative disk files

print('primer_der.eps', '-depsc');

print('primer_der.jpg', '-djpeg');

print('primer_der.tif', '-dtiff');







