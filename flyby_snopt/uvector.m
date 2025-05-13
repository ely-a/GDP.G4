function vhat = uvector(x)

% unit vector

% input

%  x = 3-component input vector

% output

%  vhat = 3-component unit vector

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vmag = sqrt(dot(x, x));

if (vmag > 0)
   vhat = x / vmag;
else
   vhat = 1.0e99;
end

