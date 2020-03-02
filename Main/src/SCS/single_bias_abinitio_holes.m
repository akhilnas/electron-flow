function single_bias_abinitio_holes(input_file)

OutGenerate=generate(input_file);

% Choice of Optimization of Fitting Factor (choice = 0) or 
% Fitting Factor Sweep (choice =1)
choice = 0;

% Choice of Width Simulation
% 0 for 8.7 nm
% 1 for 5.4 nm
% 2 for 3.2 nm
wchoice = 2;

H=getH(OutGenerate);
H=sparse(H);

if wchoice == 0
% Ab-initio data for 8.7 nm
  Extra_hole = [0 ,  0.00077778,  0.00155556,  0.00233333,  0.00311111, 0.00388889,  0.00466667,  0.00544444,  0.00622222,  0.007  ];
  Gate_Bias = [0 , -0.45404   , -0.49566699, -0.52946154, -0.55982901,-0.58817365, -0.61513207, -0.64105449, -0.66615823, -0.6905854];
%len = '8.7 nm';
elseif wchoice == 1
% Ab-initio data for 5.46 nm
Extra_hole = [0        ,  0.00020408,  0.00040816,  0.00061224,  0.00081633,...
        0.00102041,  0.00122449,  0.00142857,  0.00163265,  0.00183673,...
        0.00204082,  0.0022449 ,  0.00244898,  0.00265306,  0.00285714,...
        0.00306122,  0.00326531,  0.00346939,  0.00367347,  0.00387755,...
        0.00408163,  0.00428571,  0.0044898 ,  0.00469388,  0.00489796,...
        0.00510204,  0.00530612,  0.0055102 ,  0.00571429,  0.00591837,...
        0.00612245,  0.00632653,  0.00653061,  0.00673469,  0.00693878,...
        0.00714286,  0.00734694,  0.00755102,  0.0077551 ,  0.00795918,...
        0.00816327,  0.00836735,  0.00857143,  0.00877551,  0.00897959,...
        0.00918367,  0.00938776,  0.00959184,  0.00979592,  0.01 ];
Gate_Bias = [0        , -0.4173284 , -0.44055949, -0.45634295, -0.4690758 ,...
       -0.48013489, -0.49013216, -0.49939454, -0.50811706, -0.51642763,...
       -0.52441115, -0.5321289 , -0.53962497, -0.54693602, -0.55408766,...
       -0.56110108, -0.56799335, -0.57477842, -0.58146785, -0.58807135,...
       -0.5945971 , -0.6010521 , -0.60744236, -0.61377312, -0.62004892,...
       -0.62627374, -0.63245113, -0.63858421, -0.64467575, -0.65072826,...
       -0.65674397, -0.66272491, -0.66867288, -0.67458955, -0.6804764 ,...
       -0.68633489, -0.69216729, -0.69797131, -0.70375192, -0.70950869,...
       -0.71524216, -0.72095346, -0.72664326, -0.73231269, -0.73796407,...
       -0.74359414, -0.74920576, -0.75479927, -0.76037521, -0.76593405];
%len = '5.4 nm';
elseif wchoice == 2
% Ab-initio data for 3.2 nm
Extra_hole = [0.        ,  0.00020408,  0.00040816,  0.00061224,  0.00081633,...
        0.00102041,  0.00122449,  0.00142857,  0.00163265,  0.00183673,...
        0.00204082,  0.0022449 ,  0.00244898,  0.00265306,  0.00285714,...
        0.00306122,  0.00326531,  0.00346939,  0.00367347,  0.00387755,...
        0.00408163,  0.00428571,  0.0044898 ,  0.00469388,  0.00489796,...
        0.00510204,  0.00530612,  0.0055102 ,  0.00571429,  0.00591837,...
        0.00612245,  0.00632653,  0.00653061,  0.00673469,  0.00693878,...
        0.00714286,  0.00734694,  0.00755102,  0.0077551 ,  0.00795918,...
        0.00816327,  0.00836735,  0.00857143,  0.00877551,  0.00897959,...
        0.00918367,  0.00938776,  0.00959184,  0.00979592,  0.01 ]; 
Gate_Bias = [0.        , -0.45879268, -0.48130649, -0.49638224, -0.50841216,...
       -0.51877275, -0.52807711, -0.53665232, -0.54469359, -0.55232684,...
       -0.55963806, -0.56668856, -0.57352364, -0.58017767, -0.58667739,...
       -0.59304398, -0.59929452, -0.605443  , -0.61150097, -0.61747814,...
       -0.6233827 , -0.62922167, -0.63500105, -0.64072607, -0.64640128,...
       -0.65203068, -0.65761779, -0.66316573, -0.6686773 , -0.67415498,...
       -0.679601  , -0.68501739, -0.69040597, -0.6957684 , -0.70110618,...
       -0.70642069, -0.7117132 , -0.71698487, -0.72223674, -0.72746982,...
       -0.73268501, -0.73788314, -0.74306501, -0.74823135, -0.75338282,...
       -0.75852007, -0.7636437 , -0.76875427, -0.77385229, -0.77893826];
%len = '3.9 nm';
end
% Pre - Initialization

% Calculation of Electric Field

Disp_Electric_Field = -(Extra_hole*1.6*10^-19)/(2*3.84*3.84*10^-20) ;

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
            
             scs(H,OutGenerate,OutGenerate.surface.bias{i}); 
%             scs_new(H,OutGenerate,OutGenerate.surface.bias{i});            
            
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
    h = (x(end) - x(1))*(0.0529)*(10^-7) / length(x) ;
    Integral = 0;
    for k=1:2:(length(x)-2)
    Integral =  Integral + (h/3) * (SP.p_charge(k) + 4*SP.p_charge(k+1) + SP.p_charge(k+2)) * (1/((5.29*10^-9)^3));
    end
    SP_Extra_hole(counter1) = Integral*(3.84*3.84*10^-16) ;
    
  if(choice == 0)
    
   OutGenerate.surface_potential =  Vg(counter1);
   [Vbetterguess, A]   = DG_initial(OutGenerate);
  if wchoice == 0
  fitting_factor_initial = 6;
  elseif wchoice == 1  
  fitting_factor_initial = 2;  
  elseif wchoice == 2
  fitting_factor_initial = 0.4;
  end
  lb=0;
  ub=100;
  options.TolFun=1e-09;
  options.TolX=1e-09;
  options.MaxIter=100;
  options.MaxFunEvals=1000;
  fitting_factor_final(counter1) = lsqnonlin(@abinitio_fitting_holes,fitting_factor_initial,lb,ub,options,OutGenerate,Extra_hole(counter1),Vbetterguess,A);
  fprintf('DG fitting parameter = %f',fitting_factor_final); 
  
  % For DG Plotting purpose
    DG = open(['output/trial/dev1_sp/mat_files/file_' num2str(Vg(counter1)) '_DG.mat']);
    % Processing
    x = OutGenerate.x;
    h = (x(end) - x(1))*(0.0529)*(10^-7) / (length(x)-1) ;
    Integral = 0;
    for k=1:2:(length(x)-2)
    Integral =  Integral + (h/3) * (DG.p_charge_DG(k) + 4*DG.p_charge_DG(k+1) + DG.p_charge_DG(k+2)) * (1/((5.29*10^-9)^3));
    end
    DG_Extra_hole(counter1) = Integral*(3.84*3.84*10^-16) ;
  end
end
  
  % Plotting
  figure(1);
  plot(Gate_Bias,Extra_hole,'o');
  hold on;
  
  figure(1);
  plot(Gate_Bias,SP_Extra_hole,'--');
  hold on;
  
  if(choice == 0)
  
      figure(1);
      plot(Gate_Bias,DG_Extra_hole,'-','LineWidth',2);
      xlabel('Gate Bias(V)','FontSize', 20);
      ylabel('Extra hole(number)','FontSize', 20);
      legend('abinitio','Eff-mass S.P','D-G');
      title('Plot of Extra hole vs. Gate Bias for different methods','FontSize', 20);
%       str = {'Fin Length = ' len ,'No Doping','Oxide Sio2 of 1nm width'};
%       text(-0.3,0.007,str,'FontSize', 20);
      hold on;
      
      figure(2)
      plot(Gate_Bias,fitting_factor_final,'-','LineWidth',2);
      hold on;
      plot(Gate_Bias,fitting_factor_final,'o');
      xlabel('Gate Bias(V)','FontSize', 20);
      ylabel('DG Fitting Factor','FontSize', 20);
%       str = {'Fin Length = 8.7 nm','No Doping','Oxide Sio2 of 1nm width'};
%       text(-0.2,1.9,str,'FontSize', 20);
      title('Plot of Variation of Density-Gradient Fitting Factor with Applied Gate Voltage','FontSize', 20);
  
  elseif(choice == 1)
      fitting_factor_final = linspace(0 , 9 , 10 );
      for counter2 = 1:length(fitting_factor_final)
          for counter1 = 1:length(Extra_hole)
              % Running the Density Gradient Algorithm
              OutGenerate.surface_potential =  Vg(counter1);
              [Vbetterguess, A]   = DG_initial(OutGenerate);              
              [~] = scsdg(fitting_factor_final(counter2),OutGenerate,Vbetterguess,A);
              
              % Opening the Density Gradient File
              DG = open(['output/trial/dev1_sp/mat_files/file_' num2str(Vg(counter1)) '_DG.mat']);
              
              % Calculaing P charge integral using simpsons rule              
              x = OutGenerate.x;
              h = (x(end) - x(1))*(0.0529)*(10^-7) / (length(x)-1) ;
              Integral = 0;
              for k=1:2:(length(x)-2)
                Integral =  Integral + (h/3) * (DG.p_charge_DG(k) + 4*DG.p_charge_DG(k+1) + DG.p_charge_DG(k+2)) * (1/((5.29*10^-9)^3));
              end
              DG_Extra_hole(counter1) = Integral*(3.84*3.84*10^-16);
          end
          
          % Plotting          
          figure(1);
          plot(Gate_Bias,DG_Extra_hole,'-','LineWidth',2);
          hold on;
          
      end
      
      % Plot Specifications  
      figure(1);
      xlabel('Gate Bias(V)');
      ylabel('Extra hole (number)');
      legend('ab-initio','Effective Mass S-P','fitting factor = 0' ,'fitting factor = 1', 'fitting factor = 2', 'fitting factor = 3', 'fitting factor = 4', 'fitting factor = 5', 'fitting factor = 6', 'fitting factor = 7', 'fitting factor = 8' ,'fitting factor = 9' );
      title('Plot of Extra Hole vs. Gate Bias for different methods');
%       str = {'Fin Length = 8.7 nm','No Doping','Oxide Sio2 of 1nm width'};
%       text(-0.3,0.01,str);
  
  end
  
  
  
end