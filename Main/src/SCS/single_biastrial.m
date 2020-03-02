% Newton - Raphson implementation of CVsimulator
% Started on 21 May 2014 by dhirendra

% Maintained by : Dhirendra (dhirendra22121987@gmail.com)

function single_bias(input_file)

OutGenerate=generate(input_file);


H=getH(OutGenerate);
H=sparse(H);

% Bias Condition
OutGenerate.surface_potential = 0.05;

global Q eps0 Kb;


if ( OutGenerate.control.Solvers.SP || OutGenerate.control.Solvers.P )
    for i = 1:size(OutGenerate.surface.bias,2)
        try
            
            scs(H,OutGenerate,OutGenerate.surface.bias{i});
            Schrodinger_Poisson  = open('output/trial/dev1_sp/mat_files/file_0.05.mat');
            [ V0 ] = Poisson(Schrodinger_Poisson.n_charge,Schrodinger_Poisson.p_charge,OutGenerate);
            Ec               = OutGenerate.Ec;
            Ev               = OutGenerate.Ev;
            No_of_nodes      = length(OutGenerate.x);
            % Calculation for holes
            OutSchr_common = Schrodinger4_holes(OutGenerate,(Ec-Q*V0),(Ev-Q*V0));
            % Calculation for electrons
            mn_t      = OutGenerate.layer.layer_material.mn_t*ones(1,No_of_nodes);
            OutSchr_t = Schrodinger4_electrons(OutGenerate,(Ec-Q*V0),(Ev-Q*V0),mn_t);
            OutSchr_t.Psi2_V = OutSchr_common.Psi2_V;
            OutSchr_t.Eigen_val_V = OutSchr_common.Eigen_val_V;
            OutSchr_t.direction = 't';
            [rho_t, n_charge_t, p_charge, rho1] = charge_density(V0,OutGenerate,(Ec-Q*V0),(Ev-Q*V0),OutSchr_t,0);
            mn_l      = OutGenerate.layer.layer_material.mn_l*ones(1,No_of_nodes);
            OutSchr_l = Schrodinger4_electrons(OutGenerate,(Ec-Q*V0),(Ev-Q*V0),mn_l);
            OutSchr_l.Psi2_V = OutSchr_common.Psi2_V;
            OutSchr_l.Eigen_val_V = OutSchr_common.Eigen_val_V;
            OutSchr_l.direction = 'l';
            [rho_l, n_charge_l, p_charge, rho1] = charge_density(V0,OutGenerate,(Ec-Q*V0),(Ev-Q*V0),OutSchr_l,0);
            n_charge = n_charge_t + n_charge_l;
            figure(1);
            plot(Schrodinger_Poisson.OutGenerate.x*0.0529,Schrodinger_Poisson.n_charge/((5.29*10^-9)^3));
            hold all;
            plot(Schrodinger_Poisson.OutGenerate.x*0.0529,n_charge/((5.29*10^-9)^3));
            xlabel('Distance(nm)');
            ylabel('Electron Density(cm^-3)');
            
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

fprintf('\nSaving the results in:\n %s \n',output_dir);
%Generate_text_HFCV(output_dir,file_prefix);
%generate_text_output_HF(output_dir);


end
  [Vbetterguess , A]   = DG_initial(OutGenerate);
  Schrodinger_Poisson  = open('output/trial/dev1_sp/mat_files/file_0.1.mat');
  fitting_factor_initial = 30;
  lb=0;
  ub=100;
  options.TolFun=1e-30;
  options.TolX=0;
  options.MaxIter=200;
  options.MaxFunEvals=2000;
  fitting_factor_final = lsqnonlin(@DG_fitting,fitting_factor_initial,lb,ub,options,OutGenerate,Schrodinger_Poisson,Vbetterguess,A);
  fprintf('DG fitting parameter = %f',fitting_factor_final);
  
figure(1);
plot(Schrodinger_Poisson.OutGenerate.x*0.0529,Schrodinger_Poisson.n_charge/((5.29*10^-9)^3),'--','LineWidth',2);
ylabel('Electron Density(cm-3)');
xlabel('x(nm)');
hold all;
% figure(2);
% plot(Schrodinger_Poisson.OutGenerate.x,VS,'--','LineWidth',2);
% hold all;
 %ff = linspace(0.001,0.01,10);

 %for i = 1:10
scsdg(fitting_factor_final,OutGenerate,Vbetterguess,A);

fprintf('\n\nSimulation completed...\n');
DG = open('output/trial/dev1_sp/mat_files/file_0.1_DG.mat');
figure(1);
plot(DG.x*0.0529,DG.n_charge_DG/((5.29*10^-9)^3),'LineWidth',2);
legend('S.P','D.G');
hold all;
[ VSP ] = Poisson(Schrodinger_Poisson.n_charge,Schrodinger_Poisson.p_charge,OutGenerate);
[ VDG ] = Poisson(DG.n_charge_DG,DG.p_charge_DG,OutGenerate);
figure(3);
plot(DG.x*0.0529,VSP*19.2,'--','LineWidth',2);
ylabel('Potential(Volts)');
xlabel('x(nm)');
hold all;
figure(3);
plot(DG.x*0.0529,VDG*19.2,'LineWidth',2);
legend('S.P','D.G');
hold all;

% Error Plot
figure(2);
plot(DG.x*0.0529,(VSP-VDG)*19.2,'o');
ylabel('Difference in Potential(Volts)');
xlabel('x(nm)');
 %end
end