function f = tpbvp_objective_otb(x)

% OTB two point boundary value objective function

% input

%  x = current delta-v vector (kilometers/second)

% output

%  f = objective function (delta-v magnitude, kilometers/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global otype

% objective function (delta-v magnitude, kilometers/second)

if (otype == 1)
    
   % initial delta-v (flyby)
   
   f = norm(x);
   
else
    
   % total delta-v (rendezvous)
   
   f = norm(x(1:3)) + norm(x(4:6));
   
end

