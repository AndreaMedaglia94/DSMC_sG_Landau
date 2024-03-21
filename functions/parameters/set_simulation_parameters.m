function [p_sim, Psi] = set_simulation_parameters(p_sch, p_phys)

% number of particles
p_sim.N  = 5e7; 

% % % Legendre polynomials % % %

% order of expansion, nodes, and weights
p_sim.M             = 4;
p_sim.Nk            = p_sim.M + 1 ;
[p_sim.zk,p_sim.wk] = lgwt(p_sim.Nk,-1,1);
p_sim.nwk           = sum(p_sim.wk);

% Legendre polynomials of deg 0:p_sim.M in the nodes 1:p_sim.Nk
Psi = zeros(p_sim.M+1,p_sim.Nk); 
for k = 0:p_sim.M
    Psi(k+1,:) = legPoly(p_sim.zk,k);
end

% normalize the polynomials to have orthonormality
Norm=zeros(p_sim.M+1,1);
for k = 0:p_sim.M
    Norm(k+1)=sum(Psi(k+1,:).^2.*p_sim.wk')./p_sim.nwk;
    Psi(k+1,:)=Psi(k+1,:)./sqrt(Norm(k+1));
end


% % % time paramteres % % %

% final time of the simulation
if strcmp(p_sch.test, 'BKW')
    p_sim.tf  = 4;
elseif strcmp(p_sch.test, 'Trub')
    if strcmp(p_sch.pot, 'Maxwell')
        p_sim.tf  = 3;
    elseif strcmp(p_sch.pot, 'Coulomb')
        p_sim.tf  = 14;
    end
end

% time step
p_sim.dt  = 0.1; 

% parameter epsilon approximating the Boltzmann equation
p_sim.epsi = p_phys.rho * p_sim.dt; 

% total number of steps
p_sim.ntot   = ceil(p_sim.tf / p_sim.dt);

% time step for the observables
p_sim.t_obs    = p_sim.dt; 
p_sim.time_obs = 0:p_sim.t_obs:p_sim.tf;

% time step for the plot
p_sim.t_plt  = 1; 


% % % uncertain parameters for initial data % % %

% temperature along the x-y-z axis
if strcmp(p_sch.test, 'BKW')
    p_sim.Tx = 0.95 + 0.1 .* ((p_sim.zk+1)./2) ;
    p_sim.Ty = 0.95 + 0.1 .* ((p_sim.zk+1)./2) ;
    p_sim.Tz = 0.95 + 0.1 .* ((p_sim.zk+1)./2) ; 

    p_sim.BKW_B = 1/8;
    p_sim.BKW_C = 2/5;
    p_sim.BKW_t = 0  ;
elseif strcmp(p_sch.test, 'Trub')
    if strcmp(p_sch.pot, 'Maxwell')
        p_sim.Tx = 0.085 ;
        p_sim.Ty = 0.085 ;
        p_sim.Tz = 0.04 ; 
    elseif strcmp(p_sch.pot, 'Coulomb')
        p_sim.Tx = 0.07 ;
        p_sim.Ty = 0.07 ;
        p_sim.Tz = 0.05 ; 
    end
end

% total initial temperature
p_sim.Ttot = (p_sim.Tx + p_sim.Ty + p_sim.Tz) / 3 ;

% % % parameters for the reconstruction % % %

% right boundary of the v-domain for the reconstruction
p_sim.L     = 5;

% number of bins per direction for the histogram reconstruction
p_sim.Nbins  = 50;
p_sim.Nexact = 100;

% edges and step for the histogram reconstruction 
p_sim.VEdges = linspace(-p_sim.L, p_sim.L, p_sim.Nbins + 1);
p_sim.dV     = p_sim.VEdges(2) - p_sim.VEdges(1) ;
p_sim.VCells = p_sim.VEdges(1:end-1) + p_sim.dV/2 ;

% edges and step for the exact solution
p_sim.VEdges_exact = linspace(-p_sim.L, p_sim.L, p_sim.Nexact+1);
p_sim.dV_exact     = p_sim.VEdges_exact(2) - p_sim.VEdges_exact(1) ;
p_sim.VCells_exact = p_sim.VEdges_exact(1:end-1) + p_sim.dV_exact/2 ;
[p_sim.Vx_exact,p_sim.Vy_exact,p_sim.Vz_exact] = meshgrid(p_sim.VCells_exact,p_sim.VCells_exact,p_sim.VCells_exact);

% % % other parameters % % %

% NonLinear equation for A
if strcmp(p_sch.kernel, 'D1')
    if strcmp(p_sch.pot, 'Maxwell')
        nonlineq   = @(x) (coth(x)-1./x-exp(-p_sim.epsi)) ;
        options    = optimset('TolX',1e-24);
        p_sim.A    = fzero(nonlineq, 1./p_sim.epsi,options);
    elseif strcmp(p_sch.pot, 'Coulomb')
        p_sim.nonlineq =  @(x,y) (coth(x)-1/x-exp(-y)) ;
    end
end

end