I = diag([317500, 317500, 186000]);
rp = 6678;
ra = 48000;
e = (ra - rp) / (ra + rp);
mu = 398600;
h = sqrt(rp * mu * (1+e));
total = 0;
for theta = linspace(0, 2*pi, 100)
    r = h^2 / mu * (1 ./ (1 + e*cos(theta))) .* [e*cos(theta) sin(theta) 0]';
    t = 3 * mu ./ norm(r)^5 * norm(r)^2 * 317500; %cross(r, I*r);
    total = total + t / 100;
end