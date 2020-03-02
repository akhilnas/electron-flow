function [ Difference ] = DG_fitting(fitting_factor,OutGenerate, Schrodinger_Poisson,Vbetterguess,A)
% Function to fit the Density-Gradient solution to Schrodinger-Poisson solution
% Based on the value of fitting factor

% Calling Density-Gradient Solver
[DG_n_charge, DG_p_charge] = scsdg(fitting_factor,OutGenerate,Vbetterguess,A);

% The Paramter to be optimized i.e be equal to zero
if OutGenerate.surface_potential >=0
    Difference = norm(Schrodinger_Poisson.n_charge  - DG_n_charge )*1e06;
else
    Difference = norm(Schrodinger_Poisson.p_charge  - DG_p_charge )*1e06;
end


end

