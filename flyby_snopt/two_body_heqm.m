function ydot = two_body_heqm (~, y)

% first-order two-body heliocentric equations of motion

% input

%  y = heliocentric state vector (kilometers and kilometers/second)

% output

%  ydot = integration vector (kilometers/second/second)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global smu

rmag = norm(y(1:3));

r3 = rmag * rmag * rmag;

% acceleration due to the point-mass sun (kilometers/second/second)

acc(1) = -smu * y(1) / r3; 

acc(2) = -smu * y(2) / r3;

acc(3) = -smu * y(3) / r3;

% compute total acceleration vector (kilometers/second/second)

ydot(1) = y(4);

ydot(2) = y(5);

ydot(3) = y(6);

ydot(4) = acc(1);

ydot(5) = acc(2);

ydot(6) = acc(3);

ydot = ydot';



