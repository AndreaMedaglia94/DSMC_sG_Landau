function [obs, distr] = Solve(vx_hat, vy_hat, vz_hat, p_sch, p_sim, p_phys, obs, distr, Psi)

fprintf('Solving the equation: \n');

% time counters
counter_obs = 1 ;
counter_plt = 1 ;

% time cycle
for n=1:p_sim.ntot
    
    fprintf('t=%f\n',n*p_sim.dt);

    [vx_hat, vy_hat, vz_hat] = NanbuBabovski(vx_hat, vy_hat, vz_hat, p_sch, p_sim, p_phys, Psi);
    
    if mod(n*p_sim.dt,p_sim.t_obs)==0
        counter_obs = counter_obs + 1 ;
        [obs{1,counter_obs}, obs{2,counter_obs}, obs{3,counter_obs}, obs{4,counter_obs}, obs{5,counter_obs}]= Observables(vx_hat, vy_hat, vz_hat, p_sim, Psi);
    end
    
    if mod(n*p_sim.dt,p_sim.t_plt)==0
        if strcmp(p_sch.test, 'BKW')
            counter_plt = counter_plt + 1 ;
            [distr{1,counter_plt}, distr{2,counter_plt}, distr{3,counter_plt}, distr{4,counter_plt}] = reconstruction(vx_hat, vy_hat, vz_hat, p_sim, p_phys, Psi) ;
        end
    end

 end    
    
fprintf('Done \n');
fprintf('----------------------------------------------------------------------- \n\n\n');
end