function [n_charge_DG,p_charge_DG ] =  scsdg(fitting_factor,OutGenerate,Vbetterguess,A)
% Function that calls the Density-Gradient Solver and writes the solution data to file

global Q eps0 Kb;
Vt=Kb*OutGenerate.Tempr;

% Boundary Condition
V_gate = OutGenerate.surface_potential;


%Initialisation of Variables
output_dir       = OutGenerate.control.Output_directory;
file_prefix      = OutGenerate.control.file_prefix;
x                = OutGenerate.x;
%Making and Saving of file in directory
out_dir{1}=OutGenerate.control.Output_directory;
if output_dir(end)=='/'
    if ~exist([output_dir(1:end-1) '_HF'],'dir')
        mkdir([output_dir(1:end-1) '_HF']);
    end
    out_dir{2}=[output_dir(1:end-1) '_HF'];
else
    mkdir([output_dir '_HF']);
    out_dir{2}=[output_dir '_HF'];
end

 %DG implementation
        
        %DG Solver
         [ VDG , n_charge_DG , p_charge_DG ] = DG_Solver(fitting_factor,OutGenerate,Vbetterguess,A,Vt,Q);


%Writing Data into file
variables_name = {'VDG' , 'n_charge_DG' , 'p_charge_DG', 'x' , 'fitting_factor'};

if ~exist(fullfile(output_dir,'mat_files'),'dir')
    mkdir(fullfile(output_dir,'mat_files'))
end

filename = fullfile(output_dir,'mat_files',[file_prefix num2str(V_gate) '_DG.mat']);
save(filename,variables_name{:});
OutGenerate.control.Output_directory=out_dir{1};

end