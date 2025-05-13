% lambert4_otb.m        January 26, 2024

% NLP solution of the j2-perturbed Earth orbit Lambert problem

% Optimization Toolbox version

% flyby and rendezvous trajectory types

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

global mu neq tetol rkcoef

global otype ri vi tof rtarget vtarget drf dvf

sv1 = zeros(6, 1);

sv2 = zeros(6, 1);

rfp = zeros(3, 1);

vfp = zeros(3, 1);

% user-defined astrodynamic and utility constants

om_constants;

% initialize rkf78 algorithm

rkcoef = 1;

% rkf78 convergence tolerance

tetol = 1.0e-10;

% number of differential equations

neq = 6;

clc; home;

fprintf('\n           program lambert4 - OTB\n');

fprintf('\n< j2-perturbed Earth orbit Lambert problem >\n');

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

revmax = 0;

% compute state vectors of initial and final orbits
% (kilometers and kilometers/second)

[ri, vi] = orb2eci(mu, oev1);

[rf, vf] = orb2eci(mu, oev2);

% save the final position and velocity vectors
% (kilometers and kilometers/second)

rtarget = rf;

vtarget = vf;

% compute initial guess for delta-v vector
% using Gooding Lambert algorithm (kilometers/second)

sv1(1:3) = ri(1:3);

sv1(4:6) = vi(1:3);

sv2(1:3) = rf(1:3);

sv2(4:6) = vf(1:3);

[vito, vfto] = lambert_gooding(mu, sv1, sv2, tof, revmax);

oev3 = eci2orb1(mu, ri, vito);

% initial guess for delta-v vectors (kilometers/second)

xg(1) = vito(1) - vi(1);
xg(2) = vito(2) - vi(2);
xg(3) = vito(3) - vi(3);

if (otype == 2)
    
    % rendezvous
    
    xg(4) = vf(1) - vfto(1);
    xg(5) = vf(2) - vfto(2);
    xg(6) = vf(3) - vfto(3);
    
end

fprintf('\ntwo-body guess for initial delta-v vector and magnitude\n');

fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * xg(1));

fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * xg(2));

fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * xg(3));

fprintf('\n\ndelta-v magnitude           %12.6f  meters/second\n', 1000.0 * norm(xg(1:3)));

if (otype == 2)

    % rendezvous
    
    fprintf('\ntwo-body guess for final delta-v vector and magnitude\n');
    
    fprintf('\nx-component of delta-v      %12.6f  meters/second', 1000.0 * xg(4));
    
    fprintf('\ny-component of delta-v      %12.6f  meters/second', 1000.0 * xg(5));
    
    fprintf('\nz-component of delta-v      %12.6f  meters/second', 1000.0 * xg(6));
    
    fprintf('\n\ndelta-v magnitude           %12.6f  meters/second', 1000.0 * norm(xg(4:6)));
    
    fprintf('\n\ntotal delta-v               %12.6f  meters/second\n\n', ...
        1000.0 * (norm(xg(1:3)) + norm(xg(4:6))));
    
end

% pre-allocate lower and upper bounds vectors

if (otype == 1)
    
    % flyby

    xlwr = zeros(3, 1);
    
    xupr = zeros(3, 1);
    
else
    
    % rendezvous

    xlwr = zeros(6, 1);
    
    xupr = zeros(6, 1);
    
end

% define lower and upper bounds for components of delta-v vectors
% (kilometers/second)

for i = 1:1:3
    
    xlwr(i) = min(- 1.1 * norm(xg(1:3)), - 0.75);
    
    xupr(i) = max(+ 1.1 * norm(xg(1:3)), + 0.75);
    
end

if (otype == 2)
    
    % rendezvous
    
    for i = 4:1:6
        
        xlwr(i) = min(-1.1 * norm(xg(4:6)), - 0.75);
        
        xupr(i) = max(+1.1 * norm(xg(4:6)), + 0.75);
        
    end
    
end

% solve trajectory optimization problem

options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', ...
    'MaxFunctionEvaluations', 5000);

[x, fval] = fmincon('tpbvp_objective_otb', xg, [], [], [], [], ...
    xlwr, xupr, 'tpbvp_constraints_otb', options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% predict the final mission orbit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi(1) = ri(1);
xi(2) = ri(2);
xi(3) = ri(3);

xi(4) = vi(1) + x(1);
xi(5) = vi(2) + x(2);
xi(6) = vi(3) + x(3);

% initial guess for step size (seconds)

h = 10.0;

% initial time (seconds)

ti = 0.0;

% final time (seconds)

tf = tof;

% integrate equations of motion

xf = rkf78('j2eqm', neq, ti, tf, h, tetol, xi);

for i = 1:1:3
    
    % initial velocity vector of the transfer orbit
    
    vito(i) = vi(i) + x(i);
    
    % predicted state vector of the final mission orbit
    
    rfp(i) = xf(i);
    
    if (otype == 2)
        
        % rendezvous

        vfp(i) = xf(i + 3) + x(i + 3);
        
    else
        
        % flyby

        vfp(i) = vf(i);
        
    end
    
end

fprintf('\n\n           program lambert4 - OTB\n');

fprintf('\n< j2-perturbed Earth orbit Lambert problem >\n');

fprintf('\norbital elements and state vector of the initial orbit');
fprintf('\n------------------------------------------------------\n');

oeprint1(mu, oev1);

svprint(ri, vi);

fprintf('orbital elements and state vector of the transfer orbit after the initial delta-v');
fprintf('\n---------------------------------------------------------------------------------\n');

oev = eci2orb1 (mu, ri, vito);

oeprint1(mu, oev);

svprint(ri, vito);

fprintf('orbital elements and state vector of the transfer orbit prior to the final delta-v');
fprintf('\n----------------------------------------------------------------------------------\n');

oev = eci2orb1 (mu, rfp, (xf(4:6))');

oeprint1(mu, oev);

svprint(rfp', (xf(4:6))');

fprintf('orbital elements and state vector of the final orbit');
fprintf('\n----------------------------------------------------\n');

oev = eci2orb1 (mu, rfp, vfp);

oeprint1(mu, oev);

svprint(rfp, vfp);

fprintf('initial delta-v vector and magnitude');
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
