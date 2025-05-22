% lambert3.m     January 26, 2024

% j2 perturbed solution of the Earth orbit Lambert problem

% shooting method with state transition matrix updates

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;

global mu j2 rkcoef

om_constants;

% initialize rkf78 integrator

rkcoef = 1;

clc; home;

fprintf('\n                    program lambert3\n');

fprintf('\n        j2 perturbed Earth orbit lambert problem \n');

fprintf('\n  shooting method with state transition matrix updates \n');

% begin simulation

fprintf('\nclassical orbital elements of the initial orbit\n');

oev1 = getoe([1;1;1;1;1;1]);

fprintf('\nclassical orbital elements of the final orbit \n');

oev2 = getoe([1;1;1;1;1;1]);

% request transfer time (minutes)

while(1)
    
   fprintf('\nplease input the transfer time in minutes\n');
   
   ttmins = input('? ');
   
   if (ttmins > 0.0)

      break;

   end
end   

% time of flight (seconds)

tof = 60.0 * ttmins;

direct = 1;

revmax = 0;

% compute state vectors of initial and final orbits
% (kilometers and kilometers/second)

[ri, vi] = orb2eci(mu, oev1);
      
[rf, vf] = orb2eci(mu, oev2);

% save the final position and velocity vectors 
% (kilometers and kilometers/second)

rfsaved = rf;

vfsaved = vf;

% compute initial guess for initial delta-v vector (kilometers/second)

sv1 = zeros(6, 1);

sv2 = zeros(6, 1);

sv1(1:3) = ri(1:3);

sv1(4:6) = vi(1:3);

sv2(1:3) = rf(1:3);

sv2(4:6) = vf(1:3);

[vito, vfto] = lambert_gooding(mu, sv1, sv2, tof, revmax);

% two-body transfer orbital elements after first impulse

oev3 = eci2orb1(mu, ri, vito);

% initial guess for initial delta-v vector (kilometers/second)

dvi(1) = vito(1) - vi(1);
dvi(2) = vito(2) - vi(2);
dvi(3) = vito(3) - vi(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine j2-perturbed solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tetol = 1.0e-10;

neq = 6;

niter = 0;

vnew = vito;

while(1)
    
   niter = niter + 1;
     
   % compute state transition matrix
   
   [rf, vf, stm] = stm2(mu, tof, ri, vnew);
   
   % extract velocity sub-matrix
   
   stm12(1:3, 1:3) = stm(1:3, 4:6);
      
   % load current state vector after first impulse
   % (kilometers and kilometers/second)
   
   x(1) = ri(1);
   x(2) = ri(2);
   x(3) = ri(3);
   
   x(4) = vnew(1);
   x(5) = vnew(2);
   x(6) = vnew(3);
   
   h = 10.0;
   
   ti = 0.0;
   
   tf = tof;
   
   % integrate equations of motion
   
   xout = rkf78('j2eqm', neq, ti, tf, h, tetol, x);
   
   % calculate final position vector difference between
   % user input and prediction (kilometers)
   
   drf(1) = rfsaved(1) - xout(1);
   drf(2) = rfsaved(2) - xout(2);
   drf(3) = rfsaved(3) - xout(3);
   
   % check for convergence
   
   drss = norm(drf);
   
   if (drss < 0.0000001 || niter > 20)
   
      break;
      
   end
   
   % compute delta-v correction (kilometers/second)
   
   dvi = stm12 \ drf';
   
   % update initial velocity vector of transfer orbit 
   % after delta-v (kilometers/second)
   
   vnew = vnew + dvi;

end

% compute initial delta-v components and magnitude (kilometers/second)
   
dvi(1) = vnew(1) - vi(1);
dvi(2) = vnew(2) - vi(2);
dvi(3) = vnew(3) - vi(3);

dvi_mag = norm(dvi);

% compute final delta-v components and magnitude (kilometers/second)
   
dvf(1) = vfsaved(1) - xout(4);
dvf(2) = vfsaved(2) - xout(5);
dvf(3) = vfsaved(3) - xout(6);

dvf_mag = norm(dvi);

fprintf('\n                    program lambert3\n');

fprintf('\n        j2 perturbed Earth orbit lambert problem \n');

fprintf('\n  shooting method with state transition matrix updates \n\n');

fprintf('\norbital elements of the initial orbit');
fprintf('\n-------------------------------------\n');

oeprint1(mu, oev1);

fprintf('\norbital elements of the final orbit');
fprintf('\n-----------------------------------\n');

oeprint1(mu, oev2);

fprintf('\ntwo-body transfer orbit');
fprintf('\n-----------------------\n');

oeprint1(mu, oev3);

fprintf('\nj2 perturbed transfer orbit');
fprintf('\n---------------------------\n');

oev4 = eci2orb1(mu, ri, vnew);

oeprint1(mu, oev4);

fprintf('\ninitial delta-v vector and magnitude\n');
   
fprintf('\nx-component of delta-v      %16.6f  meters/second', 1000.0 * dvi(1));

fprintf('\ny-component of delta-v      %16.6f  meters/second', 1000.0 * dvi(2));

fprintf('\nz-component of delta-v      %16.6f  meters/second', 1000.0 * dvi(3));

fprintf('\n\ndelta-v magnitude           %16.6f  meters/second', 1000.0 * dvi_mag);

fprintf('\n\nfinal delta-v vector and magnitude\n');
   
fprintf('\nx-component of delta-v      %16.6f  meters/second', 1000.0 * dvf(1));

fprintf('\ny-component of delta-v      %16.6f  meters/second', 1000.0 * dvf(2));

fprintf('\nz-component of delta-v      %16.6f  meters/second', 1000.0 * dvf(3));

fprintf('\n\ndelta-v magnitude           %16.6f  meters/second', 1000.0 * dvf_mag);

fprintf('\n\ntransfer time               %16.6f  minutes \n', ttmins);

fprintf('\nfinal position vector error components and magnitude\n');

fprintf('\nx-component of delta-r      %16.6f  meters', 1000.0 * drf(1));

fprintf('\ny-component of delta-r      %16.6f  meters', 1000.0 * drf(2));

fprintf('\nz-component of delta-r      %16.6f  meters', 1000.0 * drf(3));

fprintf('\n\ndelta-r magnitude           %16.6f  meters\n\n', 1000.0 * norm(drf));


