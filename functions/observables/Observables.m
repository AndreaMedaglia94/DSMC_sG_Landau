function [mean, energy, temperature, temperaturetot, moment4] = Observables(vx_hat, vy_hat, vz_hat, p_sim, Psi)

% reconstruct the velocities in the physical space (nodes)
vx = vx_hat * Psi ;
vy = vy_hat * Psi ;
vz = vz_hat * Psi ;

% compute the observables
mean           = [sum(vx,1); sum(vy,1); sum(vz,1) ] ./ p_sim.N;

energy         = 0.5 .* sum( vx.^2 + vy.^2 + vz.^2) ./ p_sim.N;
 
temperature    = [sum(vx.^2,1); sum(vy.^2,1); sum(vz.^2,1) ] ./ p_sim.N;

temperaturetot = sum(temperature,1) ./3 ;

moment4        = sum( vx.^4 + vy.^4 + vz.^4, 1) ./ p_sim.N;

end