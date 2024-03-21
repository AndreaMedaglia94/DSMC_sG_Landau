function [f3D, var3D, f1D, var1D] = reconstruction(vx_hat, vy_hat, vz_hat, p_sim, p_phys, Psi)

% reconstruct the velocities in the physical space (nodes)
vx = vx_hat * Psi ;
vy = vy_hat * Psi ;
vz = vz_hat * Psi ;

% compute the expectation and variance of the distribution
fnum     = zeros(p_sim.Nbins,p_sim.Nbins,p_sim.Nbins);
fnum2    = zeros(p_sim.Nbins,p_sim.Nbins,p_sim.Nbins);

for k=1:p_sim.Nk

    fk    = histcn([vx(:,k), vy(:,k), vz(:,k)], p_sim.VEdges, p_sim.VEdges, p_sim.VEdges) ./ p_sim.N ;
    fk    = fk ./ ( sum(fk.*p_sim.dV^3,'all') ) .* p_phys.rho ;

    fnum    = fnum   + fk    .* p_sim.wk(k);
    fnum2   = fnum2  + fk.^2 .* p_sim.wk(k);

end

fnum      = fnum    ./ p_sim.nwk ;
fnum2     = fnum2   ./ p_sim.nwk ;

var3D     = fnum2   - fnum.^2 ; 
f3D       = fnum ;

% compute marginals
f1D       = sum( sum( f3D,   3 ), 2 ) .* p_sim.dV^2 ;
var1D     = sum( sum( var3D, 3 ), 2 ) .* p_sim.dV^2 ;

end