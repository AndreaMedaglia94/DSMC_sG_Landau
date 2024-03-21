function [obs, distr] = set_observables(vx_hat, vy_hat, vz_hat, p_sim, p_phys, Psi)

[obs{1,1}, obs{2,1}, obs{3,1}, obs{4,1}, obs{5,1}]= Observables(vx_hat, vy_hat, vz_hat, p_sim, Psi);

[distr{1,1}, distr{2,1}, distr{3,1}, distr{4,1}] = reconstruction(vx_hat, vy_hat, vz_hat, p_sim, p_phys, Psi) ;

end