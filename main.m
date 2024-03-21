%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MIT License
% 
% Copyright (c) 2024 Andrea Medaglia
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear
clc

addpath(genpath('functions')); 

% Set the structs containing the parameters
[p_phys, p_sim, p_sch, Psi] = set_parameters ;

% Set initial data
[vx_hat, vy_hat, vz_hat] = InitialData(p_sch, p_sim, Psi);

% Initialize observables
[obs, distr] = set_observables(vx_hat, vy_hat, vz_hat, p_sim, p_phys, Psi);

% Solve the equation
[obs, distr] = Solve(vx_hat, vy_hat, vz_hat, p_sch, p_sim, p_phys, obs, distr, Psi) ;

% Post processing 
postprocessingplots(obs,distr,p_sim,p_phys,p_sch) ;

fprintf('---------------------------------------------------------------------------------------------------------------- \n');
fprintf('Direct Simulation Monte Carlo - stochastic Galerkin for the space homogeneous Landau equation with random inputs\n\n');
fprintf('Code written by Andrea Medaglia (last update 21/03/2024) \n');
fprintf('https://doi.org/10.1016/j.jcp.2024.112845 \n');
fprintf('---------------------------------------------------------------------------------------------------------------- \n\n');
