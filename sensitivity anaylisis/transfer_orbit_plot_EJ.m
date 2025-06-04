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

%thetas
t1 = 0:0.1:theta1;
t2 = 0:0.1:theta2;

figure
plot(r_vec(1,:),r_vec(2,:),"k");
hold on
scatter(0, 0, 100, 'filled', 'MarkerFaceColor', [255,128,0]/255);
scatter(R1*cosd(theta1),R1*sind(theta1),25, 'filled', 'MarkerFaceColor', [64,64,64]/255);
scatter(R2*cosd(theta2),R2*sind(theta2),25, 'filled', 'MarkerFaceColor', [128,128,128]/255);
plot(R1*cosd(theta),R1*sind(theta),"r");
plot(R2*cosd(theta),R2*sind(theta),"b");
hold off
legend("Transfer Orbit", "Sun","Earth","Jupiter","Earth's orbit", "Jupiter's orbit");
grid on

figure
plot(r_vec_t(1,:),r_vec_t(2,:),"k",LineWidth=1);
hold on
plot(R1*cosd(theta),R1*sind(theta),"--r", LineWidth=1);
plot(R2*cosd(theta),R2*sind(theta),"-.b",LineWidth=1);
scatter(0, 0, 100, 'filled', 'MarkerFaceColor', [255,128,0]/255);
scatter(R1*cosd(theta1),R1*sind(theta1),25, 'filled', 'MarkerFaceColor', [50,50,50]/255);
scatter(R2*cosd(theta2),R2*sind(theta2),25, 'filled', 'MarkerFaceColor', [100,100,100]/255);
plot([0,R1*cosd(theta1)],[0,R1*sind(theta1)],"k");
plot([0,R2*cosd(theta2)],[0,R2*sind(theta2)],"k");
% plot(R1*0.45*cosd(t1), R1*0.45*sind(t1),"k");
% plot(R1*0.35*cosd(t2), R1*0.35*sind(t2),"k");
plot([1.2*R2*cosd(180),1.2*R2*cosd(0)],[R2*sind(180),R2*sind(0)],"--", 'Color',[170,170,170]/255,LineWidth=0.75);
hold off
%legend("Transfer Orbit", "Earth's orbit","Jupiter's orbit","Sun","Earth", "Jupiter");
legend("Transfer Orbit", "Earth's orbit","Jupiter's orbit");
grid on
xlim([-10e8 10e8]);
ylim([-10e8 10e8]);
xlabel("X (km)");
ylabel("Y (km)");

% From the plot:
%       periapsis: 149463000
%       apoapsis: 2053310000










