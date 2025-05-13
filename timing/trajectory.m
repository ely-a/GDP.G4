clear 
clc
close all

% EMEJSN

departureDate = datetime(2030, 5, 12, 12, 0, 0); 
tofs = 10:1:400;
v_inf_sqr_depart = zeros(1,length(tofs));
v_inf_sqr_arrival = zeros(1, length(tofs));

for i = 1:length(tofs)
    result = computeInterplanetaryTransfer('Earth', 'Venus', departureDate, departureDate+days((tofs(i))));
    v_inf_sqr_depart(i) = norm(result.deltaV1)^2;
    v_inf_sqr_arrival(i) = norm(result.deltaV2)^2;
end

figure
semilogy(tofs, v_inf_sqr_depart, DisplayName= 'Earth')
hold on
semilogy(tofs, v_inf_sqr_arrival, DisplayName= 'Earth')
legend()
grid on

planets = ['Earth', 'Mars', 'Earth', 'Jupiter', 'Neptune'];