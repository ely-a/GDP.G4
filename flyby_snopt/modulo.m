function y = modulo (x)

% modulo 2 pi function

% input

%  x = argument (radians)

% output

%  y = x modulo 2 pi

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pi2 = 2 * pi;

y = x - pi2 * fix(x / pi2);

if (y < 0)
   y = y + pi2;
end


