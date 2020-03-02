function [ Difference ] = DG_fittingvector( fitting_factor,OutGenerate, Schrodinger_Poisson,Vbetterguess,A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Density_Gradient = scsdg(fitting_factor,OutGenerate,Vbetterguess,A);
Difference = (norm(Schrodinger_Poisson.n_charge  - Density_Gradient ))*1e05;


end

