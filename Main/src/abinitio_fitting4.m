function [ Difference ] = abinitio_fitting4(fitting_factor,OutGenerate,abinitio_hole,Vg,Guess )
%Fitting as per ab-initio calculations (Single Parameter)

for counter = 1:length(Vg)
    OutGenerate.surface_potential = Vg(counter);
    [~] = scsdg(fitting_factor,OutGenerate,Guess(counter).Vbetterguess,Guess(counter).A);
    
    DG = open(['output/trial/dev1_sp/mat_files/file_' num2str(Vg(counter)) '_DG.mat']);
    % Using MATLAB Trapezoidal Function to calculate Extra hole from
    % Density-Gradient Calculations
    Extra_hole(counter) = trapz(DG.x, DG.p_charge_DG - ... 
    DG.n_charge_DG)*(0.0529)*(10^-7)* (1/((5.29*10^-9)^3))*(3.84*3.84*10^-16);
end

Difference = norm(Extra_hole - abinitio_hole)*100 ;
end

