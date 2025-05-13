% lambert4_snopt.m     January 26, 2024

% SNOPT solution of the j2-perturbed Earth orbit Lambert problem

% flyby and rendezvous trajectory types

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

global mu j2 req neq tetol rkcoef

global otype ri vi tof rtarget vtarget drf dvf

om_constants;

rfp = zeros(3, 1);

vfp = zeros(3, 1);

% initialize rkf78 algorithm

rkcoef = 1;

% rkf78 convergence tolerance

tetol = 1.0e-8;

% number of differential equations

neq = 6;

clc; home;

fprintf('\n          program lambert4_snopt\n');

fprintf('\n< j2-perturbed Earth orbit Lambert problem >\n\n');

while(1)

    fprintf('\ntrajectory type (1 = flyby, 2 = rendezvous)\n');

    otype = input('? ');

    if (otype == 1 || otype == 2)

        break;

    end

end

fprintf('\nclassical orbital elements of the initial orbit\n');

oev1 = getoe([1;1;1;1;1;1]);

fprintf('\nclassical orbital elements of the final orbit \n');

oev2 = getoe([1;1;1;1;1;1]);

% transfer time

while(1)

    fprintf('\nplease input the transfer time in minutes\n');

    ttmins = input('? ');

    if (ttmins > 0.0)

        break;

    end

end

% time of flight (seconds)

tof = 60.0 * ttmins;

% compute state vectors of initial and final orbits

[ri, vi] = orb2eci(mu, oev1);

[rf, vf] = orb2eci(mu, oev2);

% save the final position and velocity vectors

rtarget = rf;

vtarget = vf;

% compute initial guess for delta-v vector
% using Gooding Lambert algorithm

sv1(1:3) = ri(1:3);

sv1(4:6) = vi(1:3);

sv2(1:3) = rf(1:3);

sv2(4:6) = vf(1:3);

revmax = 0;

[vito, vfto] = lambert_gooding(mu, sv1, sv2, tof, revmax);

oev3 = eci2orb1(mu, ri, vito');

% allocate working vectors

% ncv = number of control variables

% nmc = number of mission constraints

if (otype == 1)

    % flyby

    ncv = 3;

    nmc = 3;

    xg = zeros(ncv, 1);

    xlwr = zeros(ncv, 1);

    xupr = zeros(ncv, 1);

    flwr = zeros(nmc + 1, 1);

    fupr = zeros(nmc + 1, 1);

else

    % rendezvous

    ncv = 6;

    nmc = 6;

    xg = zeros(ncv, 1);

    xlwr = zeros(ncv, 1);

    xupr = zeros(ncv, 1);

    flwr = zeros(nmc + 1, 1);

    fupr = zeros(nmc, 1);

end

% initial guess for delta-v vectors (kilometers/second)

for i = 1:1:3

    xg(i) = vito(i) - vi(i);

end

if (otype == 2)

    % rendezvous

    for i = 1:1:3

        xg(i + 3) = vf (i) - vfto(i);

    end

end

fprintf('\n\ntwo-body guess for initial delta-v vector and magnitude\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * xg(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * xg(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * xg(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second\n', 1000.0 * norm(xg(1:3)));

if (otype == 2)

    fprintf('\ntwo-body guess for final delta-v vector and magnitude\n');

    fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * xg(4));

    fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * xg(5));

    fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * xg(6));

    fprintf('\n\ndelta-v magnitude           %12.6f  meters/second', 1000.0 * norm(xg(4:6)));

    fprintf('\n\ntotal delta-v               %12.6f  meters/second\n', ...
        1000.0 * (norm(xg(1:3)) + norm(xg(4:6))));

end

% define lower and upper bounds for components of delta-v vectors
% (kilometers/second)

xlwr(1:3) = min(-1.1 * norm(xg(1:3)), - 0.75);

xupr(1:3) = max(+1.1 * norm(xg(1:3)), + 0.75);

if (otype == 2)

    % rendezvous

    xlwr(4:6) = min(-1.1 * norm(xg(4:6)), - 0.75);

    xupr(4:6) = max(+1.1 * norm(xg(4:6)), + 0.75);

end

% bounds on objective function

flwr(1) = 0.0;
fupr(1) = +Inf;

% enforce final position vector equality constraints (kilometers)

flwr(2:4) = 0.0;
fupr(2:4) = 0.0;

if (otype == 2)

    % rendezvous - enforce final velocity vector equality constraints
    % (kilometers/second)

    flwr(5:7) = 0.0;
    fupr(5:7) = 0.0;

end

% pre-allocate SNOPT arrays according to the value of otype

if (otype == 1)

    xmul = zeros(3, 1);

    xstate = zeros(3, 1);

    fmul = zeros(4, 1);

    fstate = zeros(4, 1);

else

    xmul = zeros(6, 1);

    xstate = zeros(6, 1);

    fmul = zeros(7, 1);

    fstate = zeros(7, 1);

end

% solve the orbital TPBVP using SNOPT

snscreen on;

[x, f, inform, xmul, fmul] = snopt(xg, xlwr, xupr, xmul, xstate, ...
    flwr, fupr, fmul, fstate, 'tpbvp_snopt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predict the final mission orbit using j2-perturbed equations of motion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi(1:3) = ri(1:3);

xi(4:6) = vi(1:3) + x(1:3);

% initial guess for step size (seconds)

h = 10.0;

% initial time (seconds)

ti = 0.0;

% final time (seconds)

tf = tof;

% integrate equations of motion

xf = rkf78 ('j2eqm', neq, ti, tf, h, tetol, xi);

for i = 1:1:3

    % initial velocity vector of the transfer orbit

    vito(1:3) = vi(1:3) + x(1:3);

    % predicted state vector of the final mission orbit

    rfp(1:3) = xf(1:3);

    if (otype == 2)

        for ii = 1:1:3

            vfp(ii) = xf(ii + 3) + x(ii + 3);

        end

    else

        vfp(1:3) = vf(1:3);

    end

end

fprintf('\n\n         program lambert4_snopt\n');

fprintf('\n< j2-perturbed Earth orbit Lambert problem >\n');

fprintf('\norbital elements and state vector of the initial orbit');
fprintf('\n------------------------------------------------------\n');

oeprint1(mu, oev1);

svprint(ri, vi);

fprintf('\norbital elements and state vector of the transfer orbit after the initial delta-v');
fprintf('\n---------------------------------------------------------------------------------\n');

oev = eci2orb1 (mu, ri, vito');

oeprint1(mu, oev);

svprint(ri, vito);

fprintf('\norbital elements and state vector of the transfer orbit prior to the final delta-v');
fprintf('\n----------------------------------------------------------------------------------\n');

oev = eci2orb1(mu, rfp, (xf(4:6))');

oeprint1(mu, oev);

svprint(rfp, (xf(4:6)));

fprintf('\norbital elements and state vector of the final orbit');
fprintf('\n----------------------------------------------------\n');

oev = eci2orb1 (mu, rfp, vfp);

oeprint1(mu, oev);

svprint(rfp, vfp);

fprintf('\ninitial delta-v vector and magnitude');
fprintf('\n------------------------------------\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * x(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * x(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * x(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second\n', 1000.0 * norm(x(1:3)));

if (otype == 2)

    fprintf('\nfinal delta-v vector and magnitude');
    fprintf('\n----------------------------------\n');

    fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * x(4));

    fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * x(5));

    fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * x(6));

    fprintf('\n\ndelta-v magnitude           %12.6f  meters/second\n', 1000.0 * norm(x(4:6)));

    fprintf('\ntotal delta-v               %12.6f  meters/second\n', ...
        1000.0 * (norm(x(1:3)) + norm(x(4:6))));

end

fprintf('\nfinal position vector error components and magnitude');
fprintf('\n----------------------------------------------------\n');

fprintf('\nx-component of delta-r      %12.8f  meters', 1000.0 * drf(1));

fprintf('\ny-component of delta-r      %12.8f  meters', 1000.0 * drf(2));

fprintf('\nz-component of delta-r      %12.8f  meters', 1000.0 * drf(3));

fprintf('\n\ndelta-r magnitude           %12.8f  meters\n', 1000.0 * norm(drf));

if (otype == 2)

    % rendezvous

    fprintf('\nfinal velocity vector error components and magnitude');
    fprintf('\n----------------------------------------------------\n');

    fprintf('\nx-component of delta-v      %12.8f  meters/second', 1000.0 * dvf(1));

    fprintf('\ny-component of delta-v      %12.8f  meters/second', 1000.0 * dvf(2));

    fprintf('\nz-component of delta-v      %12.8f  meters/second', 1000.0 * dvf(3));

    fprintf('\n\ndelta-v magnitude           %12.8f  meters/second\n', 1000.0 * norm(dvf));

end

fprintf('\ntransfer time               %12.6f  minutes\n\n', ttmins);


