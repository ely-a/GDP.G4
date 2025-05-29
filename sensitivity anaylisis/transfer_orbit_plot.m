clc
clear
close all

h = 6.080792635971896e+09;
mu = 1.327e11;
e = 0.864294954010606;
theta = 0:1:360;
theta_t = (360-356.4033):0.1:137.9554;

r = ((h*h/mu)./(1+e*cosd(theta)));
r_t = ((h*h/mu)./(1+e*cosd(theta_t)));

r_vec = r .* [cosd(theta);sind(theta)];
r_vec_t = r_t .* [cosd(theta_t);sind(theta_t)];

theta1 = 356.4033;
theta2 = 137.9554;
R1 = 149.6e6;       % Distance to Earth from Sun (km)
R2 = 778e6;         % Distance to Jupiter from Sun (km)

figure
plot(r_vec(1,:),r_vec(2,:),"k");
hold on
scatter(0, 0, 100, 'filled', 'MarkerFaceColor', [255,128,0]/255);
scatter(R1*cosd(theta1),R1*sind(theta1),25, 'filled', 'MarkerFaceColor', [0,0,255]/255);
scatter(R2*cosd(theta2),R2*sind(theta2),25, 'filled', 'MarkerFaceColor', [0,255,0]/255);
plot(R1*cosd(theta),R1*sind(theta),"r");
plot(R2*cosd(theta),R2*sind(theta),"b");
hold off
legend("Transfer Orbit", "Sun","Earth","Jupiter","Earth's orbit", "Jupiter's orbit");
grid on

figure
plot(r_vec_t(1,:),r_vec_t(2,:),"k",LineWidth=1);
hold on
plot(R1*cosd(theta),R1*sind(theta),"r", LineWidth=1);
plot(R2*cosd(theta),R2*sind(theta),"b",LineWidth=1);
scatter(0, 0, 100, 'filled', 'MarkerFaceColor', [255,128,0]/255);
scatter(R1*cosd(theta1),R1*sind(theta1),25, 'filled', 'MarkerFaceColor', [0,0,255]/255);
scatter(R2*cosd(theta2),R2*sind(theta2),25, 'filled', 'MarkerFaceColor', [0,255,0]/255);
hold off
legend("Transfer Orbit", "Earth's orbit","Jupiter's orbit","Sun","Earth", "Jupiter");
grid on
xlim([-8e8 12e8]);


% From the plot:
%       periapsis: 149463000
%       apoapsis: 2053310000










