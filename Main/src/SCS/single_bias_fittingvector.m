% Newton - Raphson implementation of CVsimulator
% Started on 21 May 2014 by dhirendra

% Maintained by : Dhirendra (dhirendra22121987@gmail.com)

function single_bias(input_file)

OutGenerate=generate(input_file);


H=getH(OutGenerate);
H=sparse(H);


if ( OutGenerate.control.Solvers.SP || OutGenerate.control.Solvers.P )
    for i = 1:size(OutGenerate.surface.bias,2)
        try
            
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

fprintf('\nSaving the results in:\n %s \n',output_dir);
%Generate_text_HFCV(output_dir,file_prefix);
%generate_text_output_HF(output_dir);


end
  [Vbetterguess , A]   = DG_initial(OutGenerate);
  Schrodinger_Poisson  = open('output/trial/dev1_sp/mat_files/file_0.1.mat');
  new = open('output/trial/dev1_sp/mat_files/file_0.1_DG.mat');
  fitting_factor_initial = new.fitting_factor;
  lb=-Inf;
  ub=Inf;
  options.TolFun=1e-15;
  options.TolX=1e-60;
  options.MaxIter=1000;
  options.MaxFunEvals=1000000;
  fitting_factor_final = lsqnonlin(@DG_fittingvector,fitting_factor_initial,lb,ub,options,OutGenerate,Schrodinger_Poisson,Vbetterguess,A);
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