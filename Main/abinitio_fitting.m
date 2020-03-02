function [ Difference ] = abinitio_fitting( fitting_factor ,OutGenerate,Vg,Extra_electron)
%Fitting Ab_initio data to Density Gradient Solver

for counter1 = 1:length(Extra_electron)

OutGenerate.surface_potential =  Vg(counter1);

[Vbetterguess, A]   = DG_initial(OutGenerate);

Density_Gradient = scsdg(fitting_factor,OutGenerate,Vbetterguess,A);

% Calculaing N chagrge integral using simpsons rule

x = OutGenerate.x;
h = (x(end) - x(1))*(0.0529)*(10^-7) / length(x) ;
Integral = 0;
for k=1:2:(length(x)-2)
Integral =  Integral + (h/3) * (Density_Gradient(k) + 4*Density_Gradient(k+1) + Density_Gradient(k+2)) * (1/((5.29*10^-9)^3));
end
DG_Extra_electron(counter1) = Integral*(3.84*3.84*10^-16) ;

end
Difference = norm(Extra_electron + DG_Extra_electron);
end

