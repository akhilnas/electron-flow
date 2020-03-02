function [ Difference ] = abinitio_fitting5(fitting_factor,OutGenerate,abinitio_hole,Vbetterguess,A )
%Fitting as per ab-initio calculations

[DG_n_charge, DG_p_charge] = scsdg(fitting_factor,OutGenerate,Vbetterguess,A);

% Calculaing N chagrge integral using simpsons rule

% Using MATLAB Trapezoidal Function to calculate Extra Electron from
% Density-Gradient Calculations
Extra_hole = trapz(OutGenerate.x, DG_p_charge - ... 
DG_n_charge)*(0.0529)*(10^-7)* (1/((5.29*10^-9)^3))*(3.84*3.84*10^-16); 

Difference = ((abinitio_hole) - Extra_hole)*100 ;
end


