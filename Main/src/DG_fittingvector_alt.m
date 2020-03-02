function [ Difference ] = DG_fittingvector_alt( fitting_factor,OutGenerate, Schrodinger_Poisson,Vbetterguess,A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Density_Gradient = scsdg(fitting_factor,OutGenerate,Vbetterguess,A);
Differencecalc = Schrodinger_Poisson.n_charge  - Density_Gradient ;

Tolerence = max(abs(Differencecalc));

i = 1;

for counter = 1:length(Density_Gradient)
    if abs(Differencecalc(counter))>= (Tolerence/2)
        DG_new(i) = Density_Gradient(counter);
        SP_new(i) = Schrodinger_Poisson.n_charge(counter);
        i = i+1;
    end
end

Difference = (norm(DG_new - SP_new))*1e05;
end

