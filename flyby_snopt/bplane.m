function [bmag, bdott, bdotr, theta] = bplane(r, v)   

% b-plane coordinates

% input

%  r = planet-centered position vector
%  v = planet-centered velocity vector

% output

%  bmag  = b-plane magnitude
%  bdt   = b dot t
%  bdr   = b dot r
%  theta = b-plane angle

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu

% unit z-vector

uz(1) = 0;
uz(2) = 0;
uz(3) = 1;

% square of the velocity and radial distance

v2 = dot(v, v);   
    
rmag = sqrt(dot(r, r));    
    
% angular momentum vector    

hv = cross(r, v);
    
% angular momentum magnitude and unit vector

hmag = sqrt(dot(hv, hv));

for i = 1:1:3
    hvhat(i) = hv(i) / hmag;
end
    
% reciprocal of the semi-major axis 

ai = 2 / rmag - v2 / mu;  
    
% orbital eccentricity  

ecc2 = 1 - hmag * hmag * ai / mu;  
       
ecc = sqrt(ecc2);   
    
% periapsis radial vector (xi)   

coef1 = (v2 - mu / rmag) / (mu * ecc);

coef2 = -dot(r, v) / (mu * ecc);

for i = 1:1:3
    xi(i) = coef1 * r(i) + coef2 * v(i);
end

ximag = sqrt(dot(xi, xi));

for i = 1:1:3
    xihat(i) = xi(i) / ximag;
end

% complete the orthogonal reference frame 

eta = cross(hvhat, xihat);
    
% s-vector  

coef1 = 1 / ecc;

coef2 = sqrt(1 - 1 / ecc2);

for i = 1:1:3
    svec(i) = coef1 * xihat(i) + coef2 * eta(i);    
end    

% t-vector   

tvec = cross(svec, uz);

tvecm = sqrt(dot(tvec, tvec));

for i = 1:1:3
    tvechat(i) = tvec(i) / tvecm;
end

% r-vector    

rvec = cross(svec, tvechat);
    
% unit b-vector  

coef1 = -coef1;

for i = 1:1:3
    bvec(i) = coef2 * xihat(i) + coef1 * eta(i);    
end    

% theta angle   

bdott = dot(bvec, tvechat); 

bdotr = dot(bvec, rvec);

theta = atan3(bdotr, bdott);

% b-plane magnitude

bmag = sqrt(ecc2 - 1) / abs(ai);

bdott = bmag * bdott;

bdotr = bmag * bdotr;

