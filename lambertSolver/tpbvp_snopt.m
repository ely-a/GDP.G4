function [f, g] = tpbvp_snopt(x)

% SNOPT two point boundary value objective function
% and state vector constraints

% input

%  x = current delta-v vector (kilometers/second)

% output

%  f(1) = objective function (delta-v magnitude, kilometers/second)
%  f(2) = rx equality constraint delta (kilometers)
%  f(3) = ry equality constraint delta (kilometers)
%  f(4) = rz equality constraint delta (kilometers)
%  f(5) = vx equality constraint delta (kilometers/seconds)
%  f(6) = vy equality constraint delta (kilometers/seconds)
%  f(7) = vz equality constraint delta (kilometers/seconds)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global otype neq tetol

global ri vi tof rtarget vtarget drf dvf

if (otype == 1)

   f = zeros(4, 1);

else

   f = zeros(7, 1);

end

% load current state vector of transfer orbit
  
xi(1:3) = ri(1:3);

xi(4:6) = vi(1:3) + x(1:3);

% initial guess for step size (seconds)

h = 10.0;
  
% initial time (seconds)

ti = 0.0;

% final time (seconds)

tf = tof;
   
% integrate j2 equations of motion
   
xf = rkf78 ('j2eqm', neq, ti, tf, h, tetol, xi);

% objective function (delta-v magnitude, kilometers/second)

if (otype == 1)

   % initial delta-v only (flyby)
   
   f(1) = norm(x);

else

   % total delta-v (rendezvous)
   
   f(1) = norm(x(1:3)) + norm(x(4:6));

end

% final position vector equality constraints (kilometers)
    
f(2) = rtarget(1) - xf(1);

f(3) = rtarget(2) - xf(2);

f(4) = rtarget(3) - xf(3);

if (otype == 2)

   % final velocity vector (kilometers/second)

   vf(1) = xf(4) + x(4);
   vf(2) = xf(5) + x(5);
   vf(3) = xf(6) + x(6);

end

if (otype == 2)

   % enforce final velocity vector equality constraints (kilometers/second)

   f(5) = vtarget(1) - vf(1);

   f(6) = vtarget(2) - vf(2);

   f(7) = vtarget(3) - vf(3); 

end

% save state vector deltas for print summary

for i = 1:1:3

    drf(i) = f(i + 1);
    
    if (otype == 2)

       % rendezvous
       
       dvf(i) = f(i + 4);

    end

end

% no derivatives

g = [];
