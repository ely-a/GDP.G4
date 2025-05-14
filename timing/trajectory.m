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

min_flyby = [300, 300, 200, 1e5, 5e4, 5e4, 5e4];
planet_radii = [6052, 6378, 3390, 69911, 58232, 25362, 24622];
flyby_radii = min_flyby + planet_radii;

departureDateInitial = datetime(2030, 5, 12, 12, 0, 0); 
launch_v_inf = 7;
v_inf_sqr = launch_v_inf^2;

v_leaving_list = [];
v_arriving_list = [];
tof_list = [];

for launchDelay = 0:20:20
    total_tof = 0;
    total_dv = 0;
    disp(launchDelay)
    departureDate = departureDateInitial + launchDelay;
    for p = 1:length(planets)-1
        tof = 1;
        while tof<4000
            result = computeInterplanetaryTransfer(planets{p}, planets{p+1}, departureDate, departureDate+tof, false);
            v_inf_sqr_depart_new = norm(result.deltaV1)^2;
            turn_angle = acos(dot(result.deltaV1, result.deltaV2) / norm(result.deltaV1) / norm(result.deltaV2));
            if tof >= 2
                if min(v_inf_sqr_depart_new, v_inf_sqr_depart_old) <= v_inf_sqr && v_inf_sqr <= max(v_inf_sqr_depart_new, v_inf_sqr_depart_old)
                    tof_final = interp1([v_inf_sqr_depart_new, v_inf_sqr_depart_old], [tof, tof-1], v_inf_sqr);
                    result = computeInterplanetaryTransfer(planets{p}, planets{p+1}, departureDate, departureDate+tof_final, true);
                    v_inf_sqr = norm(result.deltaV2)^2;
                    break
                end
            end
            v_inf_sqr_depart_old = v_inf_sqr_depart_new;
            tof = tof+1;
        end
        if p == 1
            v_inf_leaving = v_inf_sqr_depart_new;
        end
        if p == length(planets)-1
            v_inf_arriving = v_inf_sqr_depart_new;
        end
        total_tof = total_tof + tof_final
        % max. turn angle
        sma = 1.327e11/v_inf_sqr_depart_new;
        ecc = 1 + flyby_radii(p)/sma;
        max_turn = 2 * asin(1/ecc) * 180/pi;
        departureDate = departureDate + days(tof_final);
    end
    v_leaving_list = [v_leaving_list v_inf_leaving];
    v_arriving_list = [v_arriving_list v_inf_arriving];
    tof_list = [tof_list total_tof];
end