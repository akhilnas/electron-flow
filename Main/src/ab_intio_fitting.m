function [ Difference ] = ab_initio_fitting( fitting_factor,OutGenerate, ab_value,Vbetterguess,A )
%Fitting as per ab-initio calculations

Density_Gradient = scsdg(fitting_factor,OutGenerate,Vbetterguess,A);

% Calculaing N chagrge integral using simpsons rule

x = OutGenerate.x;
h = (x(end) - x(1))*(0.0529)*(10^-7) / length(x) ;
Integral = 0;
for k=1:2:(length(x)-2)
Integral =  Integral + (h/3) * (n_charge(k) + 4*n_charge(k+1) + n_charge(k+2)) * (1/((5.29*10^-9)^3));
end
Extra_electron = Integral*(3.84*3.84*10^-16) ;

Difference = ab_value - Extra_electron ;
end

