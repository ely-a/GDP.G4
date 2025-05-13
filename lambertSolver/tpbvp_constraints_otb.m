function [c, ceq] = tpbvp_constraints_otb(x)

% OTB two point boundary value state vector constraints

% input

%  x = current delta-v vector (kilometers/second)

% output

%  ceq(1) = x-component of position constraint delta (kilometers)
%  ceq(2) = y-component of position constraint delta (kilometers)
%  ceq(3) = z-component of position constraint delta (kilometers)
%  ceq(4) = x-component of velocity constraint delta (kilometers/second)
%  ceq(5) = y-component of velocity constraint delta (kilometers/second)
%  ceq(6) = z-component of velocity constraint delta (kilometers/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global otype neq tetol

global ri vi tof rtarget vtarget drf dvf

% load current initial state vector of transfer orbit

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

% final position vector equality constraints (kilometers)

ceq(1) = rtarget(1) - xf(1);

ceq(2) = rtarget(2) - xf(2);

ceq(3) = rtarget(3) - xf(3);

if (otype == 2)

    % final velocity vector (kilometers/second)

    vf(1) = xf(4) + x(4);
    vf(2) = xf(5) + x(5);
    vf(3) = xf(6) + x(6);

end

if (otype == 2)

    % final velocity vector constraints (kilometers/second)

    ceq(4) = vtarget(1) - vf(1);

    ceq(5) = vtarget(2) - vf(2);

    ceq(6) = vtarget(3) - vf(3);

end

% save state vector deltas for print summary

for i = 1:1:3

    % flyby

    drf(i) = ceq(i);

    if (otype == 2)

        % rendezvous

        dvf(i) = ceq(i + 3);

    end

end

% no inequality constraints

c = [];



