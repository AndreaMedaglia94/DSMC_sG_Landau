function [f_BKW_3D, var_BKW_3D, f_BKW_1D, var_BKW_1D, M4] = BKW_analytic(p_sim, time)

if length(time)>1
    K  = p_sim.Ttot .* ( 1 - p_sim.BKW_C .* exp(-4 .* p_sim.BKW_B .* time) );
    M4      = 9 .* K .* ( 2.*p_sim.Ttot - K ) ; 
    f_BKW_3D = 0; var_BKW_3D = 0;
    f_BKW_1D = 0; var_BKW_1D = 0;
else

    % compute the analytic BKW solution at a specific time given the 3D grid v

    fnum     = zeros(p_sim.Nexact,p_sim.Nexact,p_sim.Nexact);
    fnum2    = zeros(p_sim.Nexact,p_sim.Nexact,p_sim.Nexact);
    M4       = 0 ;

    for k=1:p_sim.Nk

        K  = p_sim.Ttot(k) * ( 1 - p_sim.BKW_C .* exp(-4 .* p_sim.BKW_B .* time) );
        norm = 1./(2.*pi.*K).^(3/2);
        c1 = ( 5.*K - 3.*p_sim.Ttot(k) ) ./ ( 2.*K ) ;
        c2 = ( p_sim.Ttot(k) - K ) ./ ( 2.*K.^2 ) ;

        fk = norm .* exp( - (p_sim.Vx_exact.^2+p_sim.Vy_exact.^2+p_sim.Vz_exact.^2) ./ (2.*K) ) .* ...
            ( c1 + c2 .* (p_sim.Vx_exact.^2+p_sim.Vy_exact.^2+p_sim.Vz_exact.^2) ) ;

        fnum    = fnum   + fk    .* p_sim.wk(k);
        fnum2   = fnum2  + fk.^2 .* p_sim.wk(k);

        M4 = M4 + 9 .* K .* ( 2.*p_sim.Ttot(k) - K ) .* p_sim.wk(k);

    end

    fnum      = fnum    ./ p_sim.nwk ;
    fnum2     = fnum2   ./ p_sim.nwk ;

    var_BKW_3D   = fnum2   - fnum.^2 ; 
    f_BKW_3D     = fnum ;

    % compute marginals
    f_BKW_1D     = sum( sum( f_BKW_3D,   3 ), 2 ) .* p_sim.dV_exact.^2 ;
    var_BKW_1D   = sum( sum( var_BKW_3D, 3 ), 2 ) .* p_sim.dV_exact.^2 ;
    
    M4    =  M4 ./ p_sim.nwk ;

end

end