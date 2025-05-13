function bpprint(rca, vca, bpvector)

% print b-plane coordinates

% input



% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rtd = 180.0 / pi;

fprintf ('\nb-magnitude                %14.6f  kilometers', sqrt(bpvector(7)^2 + bpvector(8)^2));

fprintf ('\nb dot r                    %14.6f  kilometers', bpvector(8));

fprintf ('\nb dot t                    %14.6f  kilometers', bpvector(7));

fprintf ('\nb-plane angle              %14.6f  degrees', rtd * bpvector(6));

fprintf ('\nv-infinity                 %14.6f  meters/second', 1000.0 * bpvector(1));

fprintf ('\nr-periapsis                %14.6f  kilometers', bpvector(5));

fprintf ('\ndecl-asymptote             %14.6f  degrees', rtd * bpvector(2));

fprintf ('\nrasc-asymptote             %14.6f  degrees\n', rtd * bpvector(3));

sfpa = dot(rca, vca) / (norm(rca) * norm(vca));

fprintf ('\nflight path angle          %14.6f  degrees\n', rtd * asin(sfpa));

