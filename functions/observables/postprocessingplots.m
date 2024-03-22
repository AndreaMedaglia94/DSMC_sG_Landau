function postprocessingplots(obs,distr,p_sim,p_phys,p_sch)

fprintf('Plotting the results: ');

if strcmp(p_sch.test, 'BKW')

    % Marginals at fixed times of f and Var
    for i=0:p_sim.t_plt:p_sim.tf
        
        f1D   = distr{3,i+1};
        var1D = distr{4,i+1};
        [~, ~, f_BKW_1D, var_BKW_1D, ~] = BKW_analytic(p_sim, i) ;
        

        figure(i*2 + 1)
        plot(p_sim.VCells_exact, f_BKW_1D, 'k-','LineWidth',1.5)
        hold on
        plot(p_sim.VCells, f1D, 'ro','LineWidth',1.2,'MarkerSize',8)
        hold off
        legend('Exact BKW', 'DSMC', 'interpreter', 'latex', 'Location','northeast','FontSize',15)
        legend boxoff
        xlabel('$v_x$','interpreter', 'latex','FontSize',15)

        titl = sprintf('Marginal $\\textrm{E}_{\\textbf{z}}[f(v,t,\\textbf{z})]$, $t=%d$', i);
        title(titl,'interpreter', 'latex','FontSize',15)

        drawnow


        figure(i*2 + 2)
        plot(p_sim.VCells_exact, var_BKW_1D, 'k-','LineWidth',1.5)
        hold on
        plot(p_sim.VCells, var1D, 'ro','LineWidth',1.2,'MarkerSize',8)
        hold off
        legend('Exact BKW', 'DSMC', 'interpreter', 'latex', 'Location','northeast','FontSize',15)
        legend boxoff
        xlabel('$v_x$','interpreter', 'latex','FontSize',15)

        titl = sprintf('Marginal $\\textrm{Var}_{\\textbf{z}}[f(v,t,\\textbf{z})]$, $t=%d$', i);
        title(titl,'interpreter', 'latex','FontSize',15)
    
        drawnow  

    end

    M4_numerics_z   = cell2array(obs,5);
    M4_numerics     = sum( M4_numerics_z .* p_sim.wk, 1) ./ p_sim.nwk ;
    
    [~, ~, ~, ~, M4z] = BKW_analytic(p_sim, p_sim.time_obs) ;
    M4                = sum( M4z .* p_sim.wk, 1) ./ p_sim.nwk ;

    % 4th order moment
    figure(2)
    plot(p_sim.time_obs,M4_numerics,'ro','LineWidth',1.2,'MarkerSize',8)
    hold on
    plot(p_sim.time_obs,M4,'k-','LineWidth',1.5)
    hold off
    xlabel('$t$','Interpreter','latex','FontSize',15)
    legend('$\mathrm{E}_{\textbf{z}}[\mathrm{M}4(t,\textbf{z})]$ DSMC', '$\mathrm{E}_{\textbf{z}}[\mathrm{M}4(t,\textbf{z})]$ Exact','Interpreter','latex','FontSize',15,'location','best')
    legend boxoff
    title('Fourth order moment','Interpreter','latex','FontSize',15)
    ylim([(min(M4)-0.01*min(M4)) (max(M4)+0.01*max(M4))])

elseif strcmp(p_sch.test, 'Trub')

end

fprintf('Done \n \n \n');

end