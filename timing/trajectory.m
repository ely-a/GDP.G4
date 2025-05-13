clear 
clc
close all

% EMEJSN

% tofs = 10:1:400;
% v_inf_sqr_depart = zeros(1,length(tofs));
% v_inf_sqr_arrival = zeros(1, length(tofs));
% 
% for i = 1:length(tofs)
%     result = computeInterplanetaryTransfer('Earth', 'Venus', departureDate, departureDate+days((tofs(i))));
%     v_inf_sqr_depart(i) = norm(result.deltaV1)^2;
%     v_inf_sqr_arrival(i) = norm(result.deltaV2)^2;
% end

% figure
% semilogy(tofs, v_inf_sqr_depart, DisplayName= 'Earth')
% hold on
% semilogy(tofs, v_inf_sqr_arrival, DisplayName= 'Earth')
% legend()
% grid on

planets = {'Earth', 'Mars', 'Earth', 'Jupiter', 'Neptune'};
departureDate = datetime(2030, 5, 12, 12, 0, 0); 
launch_v_inf = 7;
v_inf_sqr = launch_v_inf^2;

for p = 1:length(planets)-1
    tof = 1;
    while tof<4000
        result = computeInterplanetaryTransfer(planets{p}, planets{p+1}, departureDate, departureDate+days(tof), false);
        v_inf_sqr_depart_new = norm(result.deltaV1)^2;
        if tof >= 2
            if min(v_inf_sqr_depart_new, v_inf_sqr_depart_old) <= v_inf_sqr && v_inf_sqr <= max(v_inf_sqr_depart_new, v_inf_sqr_depart_old)
                tof_final = interp1([v_inf_sqr_depart_new, v_inf_sqr_depart_old], [tof, tof-1], v_inf_sqr);
                result = computeInterplanetaryTransfer(planets{p}, planets{p+1}, departureDate, departureDate+days(tof_final), true);
                v_inf_sqr = norm(result.deltaV2)^2;
                break
            end
        end
        v_inf_sqr_depart_old = v_inf_sqr_depart_new;
        tof = tof+1;
    end
    sqrt(v_inf_sqr)
    tof_final
    departureDate = departureDate + days(tof_final)
end