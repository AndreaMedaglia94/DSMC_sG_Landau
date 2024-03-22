function [vx_hat, vy_hat, vz_hat] = InitialDataBKW(p_sim, Psi)

% sample initial data from the BKW solution by inverting the cumulative

% set some parameter to define the cumulative
K    = 1 - p_sim.BKW_C .* exp(-4 .* p_sim.BKW_B .* p_sim.BKW_t) ;
norm = 1./(2.*pi.*K).^(3/2);
c2   = ( 1 - K ) ./ ( 2.*K.^2 ) ;
b    = 2 * K ;

% sample random numbers 
xi = rand(p_sim.N,1);
% and evaluate the inverse of the cumulative on this numbers
xx     = linspace(0,6,5000);
yy     = 4.*pi.*norm.*c2.*( 3./8.*sqrt(pi).*b.^(5./2).*erf(xx./sqrt(b)) - b./4.*exp(-xx.^2./b).*(3.*b.*xx + 2.*xx.^3) );
rr     = interp1(yy,xx,xi);

% sample random angles and rotate particles
theta = acos( 1 - 2 .* rand(p_sim.N,1) ) ;
phi   = rand(p_sim.N,1) .* 2.*pi;

vxr    = rr.*sin(theta).*cos(phi); 
vyr    = rr.*sin(theta).*sin(phi); 
vzr    = rr.*cos(theta); 

% rescale the particles with the total temperature
% it works only if Tx=Ty=Tz
vx = sqrt(p_sim.Ttot') .* vxr;
vy = sqrt(p_sim.Ttot') .* vyr;
vz = sqrt(p_sim.Ttot') .* vzr;

% project the velocities in the spaces of the polynomials
vx_hat = zeros(p_sim.N,p_sim.M);
vy_hat = zeros(p_sim.N,p_sim.M);
vz_hat = zeros(p_sim.N,p_sim.M);
for h=1:p_sim.M+1
    vx_hat(:,h) = sum( vx .* Psi(h,:) .* p_sim.wk', 2 ) ./ p_sim.nwk ;
    vy_hat(:,h) = sum( vy .* Psi(h,:) .* p_sim.wk', 2 ) ./ p_sim.nwk ;
    vz_hat(:,h) = sum( vz .* Psi(h,:) .* p_sim.wk', 2 ) ./ p_sim.nwk ;
end

if sum( isnan(vx_hat),"all" )  ~= 0
    fprintf('Error in the initial sampling! \n');
    stop
elseif sum( isnan(vy_hat),"all" )  ~= 0
    fprintf('Error in the initial sampling! \n');
    stop
elseif sum( isnan(vz_hat),"all" )  ~= 0
    fprintf('Error in the initial sampling! \n');
    stop
end

end