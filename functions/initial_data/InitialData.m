function [vx_hat, vy_hat, vz_hat] = InitialData(p_sch, p_sim, Psi)

fprintf('Initializing the particles: ');

if strcmp(p_sch.test, 'BKW')
    [vx_hat, vy_hat, vz_hat] = InitialDataBKW(p_sim, Psi);
elseif strcmp(p_sch.test, 'Trub')
    % [vx_hat, vy_hat, vz_hat] = InitialDataTrubnikov(p_sim);
end

fprintf('Done \n\n');

end