% Centroid and sheet charge density calculations

% Density_Gradient Data file
DG = open('output/trial/dev1_sp/mat_files/file_0.09_DG.mat');
% Schrodinger-Poisson Data file
SP = open('output/trial/dev1_sp/mat_files/file_0.09.mat');

% Number of Nodes Calculation
No_of_nodes = length(DG.x);
% Only measured for one half of device
No_of_nodes = (No_of_nodes+1)/2;
% Pre-allocation
x = zeros(1,No_of_nodes);
n_charge_DG = zeros(1,No_of_nodes);
n_charge_SP = zeros(1,No_of_nodes);
New_DG = zeros(1,No_of_nodes);
New_SP = zeros(1,No_of_nodes);
% Assignment 
for counter1 = 1:No_of_nodes
    x(counter1)        = DG.x(counter1)*0.0529;
    n_charge_DG(counter1) = DG.n_charge_DG(counter1)/((5.29*10^-9)^3);
    n_charge_SP(counter1) = SP.n_charge(counter1)/((5.29*10^-9)^3);
end
% Centroid Calculation SP
for counter1 = 1:No_of_nodes
    New_SP(counter1) = x(counter1)*n_charge_SP(counter1);
end
Numerator_SP = trapz(New_SP,x);
Denominator_SP = trapz(n_charge_SP,x);
Centroid_SP(4) = Numerator_SP/Denominator_SP;
% Centroid Calculation DG
for counter1 = 1:No_of_nodes
    New_DG(counter1) = x(counter1)*n_charge_DG(counter1);
end
Numerator_DG = trapz(New_DG,x);
Denominator_DG = trapz(n_charge_DG,x);
Centroid_DG(4) = Numerator_DG/Denominator_DG;
% Sheet Charge Density Calculation DG
Ninv_SP(4) = trapz(n_charge_SP,x);
% Sheet Charge Density Calculation DG
Ninv_DG(4) = trapz(n_charge_DG,x);

% Error Comparison with Self-Consistent Poisson Solver
% Running Simulation
OutGenerate=generate('stack_files/trial/dev2.txt');
[Vbetterguess , A]   = DG_initial(OutGenerate);
% Maximum Error Difference between this and Schrodinger-Poisson
Vmax_Poisson = max(Vbetterguess-SP.V0)*19.2;
Vmax_DG      = max(DG.VDG - SP.V0)*19.2;
% Norm of differences
Poisson_norm = norm((Vbetterguess-SP.V0)*19.2);
DG_norm      = norm((DG.VDG - SP.V0)*19.2);

