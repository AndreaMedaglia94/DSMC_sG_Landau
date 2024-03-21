function [wwx_hat, wwy_hat, wwz_hat]=NanbuBabovski(wx_hat, wy_hat, wz_hat, p_sch, p_sim, p_phys, Psi)

% Random coupling of the particles
j  = randperm(p_sim.N); 
j1 = j(1:p_sim.N/2); 
j2 = j(p_sim.N/2+1:p_sim.N);

% Pre-interaction opinions on the projections
wxi_hat = wx_hat(j1,:); wxj_hat=wx_hat(j2,:);
wyi_hat = wy_hat(j1,:); wyj_hat=wy_hat(j2,:);
wzi_hat = wz_hat(j1,:); wzj_hat=wz_hat(j2,:);

qx_hat  = wxi_hat - wxj_hat; 
qy_hat  = wyi_hat - wyj_hat; 
qz_hat  = wzi_hat - wzj_hat; 

% Reconstruction of the velocities
wxi = wxi_hat * Psi;   wxj = wxj_hat * Psi; 
wyi = wyi_hat * Psi;   wyj = wyj_hat * Psi;
wzi = wzi_hat * Psi;   wzj = wzj_hat * Psi;

% Compute cos theta
qx = wxi - wxj; 
qy = wyi - wyj; 
qz = wzi - wzj;
qnorm = sqrt(qx.^2 + qy.^2 + qz.^2) ;

if strcmp( p_sch.pot, 'Maxwell') 
    one_over_tau = ones(p_sim.N/2,1) ;
elseif strcmp( p_sch.pot, 'Coulomb') 
    one_over_tau = 4.*pi.*(p_phys.e.^2./(4.*pi.*p_phys.m./2.*p_phys.epsi0)).^2.*p_phys.rho.*p_phys.Lambda./(qnorm.^3) ;
end


if strcmp(p_sch.kernel, 'D1')
    % if strcmp( p_sch.pot, 'Maxwell') 
    %     U   = rand(p_sim.N/2,1); 
    %     ct  = 1 ./ p_sim.A .* log( exp(-p_sim.A) + 2 .* U .* sinh(p_sim.A) );
    % elseif strcmp( p_sch.pot, 'Coulomb') 
    %     A   = getA(p_sim, one_over_tau, p_phys);
    %     U   = rand(p_sim.N/2,1);  
    %     ct  = getCost(U,A);
    % end
elseif strcmp(p_sch.kernel, 'D2')
    tau0    = p_sim.epsi ./ ( 2*p_phys.rho ) .* one_over_tau;
    nu_tau0 = (1 - 2 .* tau0) ;
    ct    = nu_tau0 .*( tau0<=1 ) + (-1) .* ( tau0>1 );
elseif strcmp(p_sch.kernel, 'D3')
    tau0  = p_sim.epsi ./ ( 2*p_phys.rho ) .* one_over_tau;
    ct    =  1 - 2.*tanh( tau0 ) ;
end

% Compute the collisions
st = sin(acos(ct));

qperp = sqrt(qy.^2 + qz.^2);
ep    = 2.*pi.*rand(p_sim.N/2,1);

hx = qperp .* cos(ep);
hy = -(qy.*qx.*cos(ep)+qnorm.*qz.*sin(ep))./qperp;
hz = -(qz.*qx.*cos(ep)-qnorm.*qy.*sin(ep))./qperp;

% collision matrices
Wx_hat = (hx.*st.*p_sim.wk') * Psi' ./ p_sim.nwk;
Wy_hat = (hy.*st.*p_sim.wk') * Psi' ./ p_sim.nwk;
Wz_hat = (hz.*st.*p_sim.wk') * Psi' ./ p_sim.nwk;

Vx_hat = (qx.*ct.*p_sim.wk') * Psi' ./ p_sim.nwk;
Vy_hat = (qy.*ct.*p_sim.wk') * Psi' ./ p_sim.nwk;
Vz_hat = (qz.*ct.*p_sim.wk') * Psi' ./ p_sim.nwk;

% compute the post-collisional velocities
wwxi_hat = wxi_hat - 0.5 .* (qx_hat - Vx_hat + Wx_hat);
wwyi_hat = wyi_hat - 0.5 .* (qy_hat - Vy_hat + Wy_hat);
wwzi_hat = wzi_hat - 0.5 .* (qz_hat - Vz_hat + Wz_hat);

wwxj_hat = wxj_hat + 0.5 .* (qx_hat - Vx_hat + Wx_hat);
wwyj_hat = wyj_hat + 0.5 .* (qy_hat - Vy_hat + Wy_hat);
wwzj_hat = wzj_hat + 0.5 .* (qz_hat - Vz_hat + Wz_hat);

wwx_hat = [wwxi_hat; wwxj_hat];
wwy_hat = [wwyi_hat; wwyj_hat];
wwz_hat = [wwzi_hat; wwzj_hat];

if sum( isnan(wwx_hat),"all" )  ~= 0
    test = 1 ;
elseif sum( isnan(wwy_hat),"all" )  ~= 0
    test = 1 ;
elseif sum( isnan(wwz_hat),"all" )  ~= 0
    test = 1 ;
end

end
