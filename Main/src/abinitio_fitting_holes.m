function [ Difference ] = abinitio_fitting_holes(fitting_factor,OutGenerate,abinitio_hole,Vbetterguess,A )
%Fitting as per ab-initio calculations
Density_Gradient = scsdg(fitting_factor,OutGenerate,Vbetterguess,A);
% Opening File
Vg = OutGenerate.surface_potential;
DG = open(['output/trial/dev1_sp/mat_files/file_' num2str(Vg) '_DG.mat']);

% Calculaing N chagrge integral using simpsons rule
x = OutGenerate.x;
h = (x(end) - x(1))*(0.0529)*(10^-7) / (length(x)-1) ;
Integral = 0;
for k=1:2:(length(x)-2)
Integral =  Integral + (h/3) * (DG.p_charge_DG(k) + 4*DG.p_charge_DG(k+1) + DG.p_charge_DG(k+2)) * (1/((5.29*10^-9)^3));
end
Extra_hole = Integral*(3.84*3.84*10^-16) ;
Difference = ((abinitio_hole) - Extra_hole) ;
end
