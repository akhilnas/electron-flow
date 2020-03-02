% Density-Gradient Solver 
% Fitting Based on Schrodinger-Poisson solver
% Started on 12 Jan 2017 by Akhil

% Maintained by : Akhil (eakhil711@gmail.com)

function single_bias(input_file)

% Reading of Parameters from txt file
OutGenerate=generate(input_file);

% Change of Effective Mass for holes
%OutGenerate.mp_eff(1,:) = 0.81;

H=getH(OutGenerate);
H=sparse(H);

% Input the Surface Bias on the Silicon substrate from user
Bias = input('Enter Surface Bias on silicon \n');

% Bias Condition - To be specified
OutGenerate.surface_potential =  Bias;

if ( OutGenerate.control.Solvers.SP || OutGenerate.control.Solvers.P )
    for i = 1:size(OutGenerate.surface.bias,2)
        try            
            % Effective Mass Schrodinger-Poisson Solver
            scs(H,OutGenerate,OutGenerate.surface.bias{i});             
        catch err
            idSegLast = regexp(err.identifier, '(?<=:)\w+$', 'match');
            if  strcmp(idSegLast,'ConvergenceFailed')
                fprintf('\n'); fprintf(err.message); fprintf('\n\n');
                return
            else
                rethrow(err);
            end
        end
    end

output_dir  = OutGenerate.control.Output_directory; 
file_prefix = OutGenerate.control.file_prefix;

fprintf('\nSaving the results in:\n %s %s \n',output_dir,file_prefix);
end

% Input the Schrodinger-Poisson solution
SP = open(['output/trial/dev1_sp/mat_files/file_' num2str(OutGenerate.surface_potential) '.mat']);

% Schrodinger-Poisson Plots
% Plot of Charge Density
figure(1);
if Bias >= 0
    plot(SP.OutGenerate.x*0.0529,SP.n_charge/((5.29*10^-9)^3),'--','LineWidth',2);
    ylabel('Electron Density(cm-3)');
    xlabel('x(nm)');
else
    plot(SP.OutGenerate.x*0.0529,SP.p_charge/((5.29*10^-9)^3),'--','LineWidth',2);
    ylabel('Hole Density(cm-3)');
    xlabel('x(nm)');  
end
hold all;

% Density-Gradient Implementation

% Initial Guess for Newton-Raphson method
[Vbetterguess , A]   = DG_initial(OutGenerate);

% Density-Gradient Optimization Scheme based on Schrodinger-Poisson solution
fitting_factor_initial = 5;
lb=0;
ub=100;
options.TolFun=1e-18;
options.TolX=1e-18;
options.MaxIter=500;
options.MaxFunEvals=2000;
fitting_factor_final = lsqnonlin(@DG_fitting,fitting_factor_initial,lb,ub,options,OutGenerate,SP,Vbetterguess,A);
fprintf('DG fitting parameter = %f',fitting_factor_final);   

% Calling Density-Gradient SOlver
[~,~] = scsdg(fitting_factor_final,OutGenerate,Vbetterguess,A);

DG = open(['output/trial/dev1_sp/mat_files/file_' num2str(OutGenerate.surface_potential) '_DG.mat']);

% Density-Gradient Plots
% Plot of Charge Density
figure(1);
if Bias >= 0
    plot(DG.x*0.0529,DG.n_charge_DG/((5.29*10^-9)^3),'LineWidth',2);
else
    plot(DG.x*0.0529,DG.p_charge_DG/((5.29*10^-9)^3),'LineWidth',2);
end
legend('Schrodinger-Poisson','Density-Gradient');
hold all;

fprintf('\n\nSimulation completed...\n');
end