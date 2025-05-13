function f = pv_dotm(ri, vi, x)

% magnitude of primer derivative vector

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu pvi pvdi

% compute state transition matrix

[~, ~, stm] = stm2(mu, x, ri, vi);

% compute magnitude of primer derivative vector

ppdot = stm * [pvi'; pvdi];

pdot = ppdot(4:6);

f = norm(pdot);
