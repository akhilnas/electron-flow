function single_bias_abiitio_single(input_file)

OutGenerate=generate(input_file);

% Choice of Optimization of Fitting Factor (choice = 0) or 
% Fitting Factor Sweep (choice =1)
choice = 0;

% Choice of Width Simulation
% 0 for 8.7 nm
% 1 for 5.4 nm
% 2 for 3.9 nm
wchoice = 0;

H=getH(OutGenerate);
H=sparse(H);

if wchoice == 0
% Ab-initio data for 8.7 nm fin width
Extra_electron = [0 -0.001, -0.002, -0.003, -0.004, -0.005, -0.006, -0.007, -0.008, -0.009, -0.01 , -0.011];
Gate_Bias = [0 0.48680108,  0.53434647,  0.57372595,  0.60941916, 0.64288517,  0.67480662,  0.70557336,  0.73543214,  0.76455412, 0.79306287,  0.82105225 ];
elseif wchoice == 1
% Data for 5.46 nm device fin width
Extra_electron = [-0 , -0.00020408, -0.00040816, -0.00061224, -0.00081633,... 
       -0.00102041, -0.00122449, -0.00142857, -0.00163265, -0.00183673,...
       -0.00204082, -0.0022449 , -0.00244898, -0.00265306, -0.00285714,...
       -0.00306122, -0.00326531, -0.00346939, -0.00367347, -0.00387755,...
       -0.00408163, -0.00428571, -0.0044898 , -0.00469388, -0.00489796,...
       -0.00510204, -0.00530612, -0.0055102 , -0.00571429, -0.00591837,...
       -0.00612245, -0.00632653, -0.00653061, -0.00673469, -0.00693878,...
       -0.00714286, -0.00734694, -0.00755102, -0.0077551 , -0.00795918,...
       -0.00816327, -0.00836735, -0.00857143, -0.00877551, -0.00897959,...
       -0.00918367, -0.00938776, -0.00959184, -0.00979592, -0.01      ];
Gate_Bias = [0.        ,  0.44933436,  0.47257258,  0.4883634 ,  0.5010991 ,...
        0.5121562 ,  0.52214786,  0.53140517,  0.54011577,  0.54840932,...
        0.5563719 ,  0.5640628 ,  0.57153152,  0.57881053,  0.5859266 ,...
        0.59290101,  0.59975093,  0.6064904 ,  0.6131311 ,  0.6196828 ,...
        0.62615379,  0.63255117,  0.63888159,  0.64514928,  0.65135941,...
        0.65751605,  0.66362283,  0.66968298,  0.67569935,  0.68167453,...
        0.68761085,  0.69351043,  0.69937516,  0.7052068 ,  0.71100695,...
        0.71677663,  0.72251815,  0.72823211,  0.73391979,  0.73958226,...
        0.7452205 ,  0.75083545,  0.75642763,  0.76199875,  0.76754836,...
        0.77307821,  0.77858859,  0.78408011,  0.78955337,  0.79500894];
elseif wchoice == 2    
% Data for 3.88 nm device fin width
Extra_electron = [-0        , -0.00020408, -0.00040816, -0.00061224, -0.00081633,...
       -0.00102041, -0.00122449, -0.00142857, -0.00163265, -0.00183673,...
       -0.00204082, -0.0022449 , -0.00244898, -0.00265306, -0.00285714,...
       -0.00306122, -0.00326531, -0.00346939, -0.00367347, -0.00387755,...
       -0.00408163, -0.00428571, -0.0044898 , -0.00469388, -0.00489796,...
       -0.00510204, -0.00530612, -0.0055102 , -0.00571429, -0.00591837,...
       -0.00612245, -0.00632653, -0.00653061, -0.00673469, -0.00693878,...
       -0.00714286, -0.00734694, -0.00755102, -0.0077551 , -0.00795918,...
       -0.00816327, -0.00836735, -0.00857143, -0.00877551, -0.00897959,...
       -0.00918367, -0.00938776, -0.00959184, -0.00979592, -0.01];
Gate_Bias = [0.        ,  0.48374189,  0.50645215,  0.52172187,  0.53394579,...
        0.54450006,  0.55399781,  0.5627661 ,  0.57100015,  0.57882589,...
        0.5863293 ,  0.5935745 ,  0.60060139,  0.60744706,  0.61413807,...
        0.62069256,  0.6271357 ,  0.63347464,  0.6397224 ,  0.64588923,...
        0.65198315,  0.65801115,  0.66397925,  0.66989268,  0.67575598,...
        0.68157314,  0.68734768,  0.69308273,  0.69878106,  0.70444517,...
        0.71007729,  0.71567942,  0.7212534 ,  0.72680087,  0.73232334,...
        0.73782219,  0.74329866,  0.74875391,  0.754189  ,  0.75960491,...
        0.76500254,  0.77038271,  0.77574622,  0.78109377,  0.78642604,...
        0.79174363,  0.79704715,  0.80233713,  0.80761407,  0.81287846];
end
% Pre - Initialization


% Calculation of Electric Field

Disp_Electric_Field = -(Extra_electron*1.6*10^-19)/(2*3.84*3.84*10^-20) ;

Electric_Field_oxide = Disp_Electric_Field/(3.9 *8.854 * 10^-12);
Potential_Drop       = Electric_Field_oxide*10^-9;
Surface_Bias         = Gate_Bias - Potential_Drop;

Vg = Surface_Bias;


for counter1 = 1:length(Vg)
% Bias Condition
OutGenerate.surface_potential =  Vg(counter1);


if ( OutGenerate.control.Solvers.SP || OutGenerate.control.Solvers.P )
    for i = 1:size(OutGenerate.surface.bias,2)
        try
            
%             scs(H,OutGenerate,OutGenerate.surface.bias{i}); 
             scs_new(H,OutGenerate,OutGenerate.surface.bias{i});            
            
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
% %Generate_text_HFCV(output_dir,file_prefix);
% %generate_text_output_HF(output_dir);

end
end

% Processing the Schrodinger Poisson Files        

  for counter1 = 1:length(Vg)   
    
    SP = open(['output/trial/dev1_sp/mat_files/file_' num2str(Vg(counter1)) '.mat']);
    % Processing
    x = OutGenerate.x;
    h = (x(end) - x(1))*(0.0529)*(10^-7) / (length(x)-1) ;
    Integral = 0;
    for k=1:2:(length(x)-2)
    Integral =  Integral + (h/3) * (SP.n_charge(k) + 4*SP.n_charge(k+1) + SP.n_charge(k+2)) * (1/((5.29*10^-9)^3));
    end
    SP_Extra_electron(counter1) = -Integral*(3.84*3.84*10^-16) ;    

  end    
  
  if wchoice == 0
  fitting_factor_initial = 1;
  elseif wchoice == 1  
  fitting_factor_initial = 1;  
  elseif wchoice == 2
  fitting_factor_initial = 1;
  end
  lb=0;
  ub=100;
  options.TolFun=1e-10;
  options.TolX=1e-25;
  options.MaxIter=500;
  options.MaxFunEvals=2000;
  fitting_factor_final(counter1) = lsqnonlin(@abinitio_fitting3,fitting_factor_initial,lb,ub,options,OutGenerate,Extra_electron,Vg);
  fprintf('DG overall fitting parameter = %f',fitting_factor_final); 
  
  for counter1 = 1:length(Vg)
  % For DG Plotting purpose
    DG = open(['output/trial/dev1_sp/mat_files/file_' num2str(Vg(counter1)) '_DG.mat']);
    % Processing
    x = OutGenerate.x;
    h = (x(end) - x(1))*(0.0529)*(10^-7) / (length(x)-1) ;
    Integral = 0;
    for k=1:2:(length(x)-2)
    Integral =  Integral + (h/3) * (DG.n_charge_DG(k) + 4*DG.n_charge_DG(k+1) + DG.n_charge_DG(k+2)) * (1/((5.29*10^-9)^3));
    end
    DG_Extra_electron(counter1) = -Integral*(3.84*3.84*10^-16) ;
  end

  % Plotting   
  figure(1);
  plot(Gate_Bias,Extra_electron,'o');
  hold on;
  
  figure(1);
  plot(Gate_Bias,SP_Extra_electron,'--');
  hold on;
      
      
  
      figure(1);
      plotyy(Gate_Bias,DG_Extra_electron,'-','LineWidth',2,Gate_Bias,Percentage_Error,'-','LineWidth',2);
      xlabel('Gate Bias(V)','FontSize', 20);
      ylabel('Extra electron (number)','FontSize', 20);
      legend('abinitio','2-mass S.P','D-G');
      title('Plot of Extra electron vs. Gate Bias for different methods','FontSize', 20);      
      hold on;
      
      figure(2)
      plot(Gate_Bias,fitting_factor_final,'-','LineWidth',2);
      hold on;
      plot(Gate_Bias,fitting_factor_final,'o');
      xlabel('Gate Bias(V)','FontSize', 20);
      ylabel('DG Fitting Factor','FontSize', 20);
      title('Plot of Variation of Density-Gradient Fitting Factor with Applied Gate Voltage','FontSize', 20);         
         
          
      
      
      % Plot Specifications  
      figure(1);
      xlabel('Gate Bias(V)');
      ylabel('Extra electron (number)');
      legend('ab-initio','2-mass S-P','fitting factor = 0','fitting factor = 1', 'fitting factor = 2', 'fitting factor = 3', 'fitting factor = 4', 'fitting factor = 5', 'fitting factor = 6', 'fitting factor = 7', 'fitting factor = 8' ,'fitting factor = 9' ,'Location','SouthWest'    );
      title('Plot of Extra Electron vs. Gate Bias for different methods');

      % Error Plot
      figure(3);
      plot(Gate_Bias,DG_Extra_electron - Extra_electron)
      xlabel('Gate Bias(V)');
      ylabel('Error Extra electron (number)');
  
  
  
end