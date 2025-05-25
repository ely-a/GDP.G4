function E = solve_kepler(M, e, tol, max_iter)
    % Inputs:
    % M        : Mean anomaly in radians
    % e        : Eccentricity (0 < e < 1)
    % tol      : Convergence tolerance (e.g., 1e-8)
    % max_iter : Maximum iterations (e.g., 100)

    % Initial guess (good for small e)
    E = M;  
    
    for k = 1:max_iter
        E_next = E - (E - e*sin(E) - M)/(1 - e*cos(E));  % fixed-point update
        if abs(E_next - E) < tol
            break;
        end
        E = E_next;
    end

    if k == max_iter
        error('Kepler solver did not converge');
    end
end