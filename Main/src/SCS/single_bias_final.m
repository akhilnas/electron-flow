% Newton - Raphson implementation of CVsimulator
% Started on 21 May 2014 by dhirendra

% Maintained by : Dhirendra (dhirendra22121987@gmail.com)

function single_bias(input_file)

OutGenerate=generate(input_file);


H=getH(OutGenerate);
H=sparse(H);

% Ab-initio data

s = load('Dg_Dhirendra/potentials.txt','-ascii');

% Calculation of Electric Field

for counter = 1:11

Electric_Field(counter) = (s(300,counter+2) - s(200,counter+2))/(s(300,1)-s(200,1)) ;
Surface_potential(counter) = Electric_Field(counter)*1 ;

end

Vg = Surface_potential ;

% Potential in Silicon

% Discarding Oxide
s(1:354,:) = [] ;
s(671:end,:) = [];

% Correction
for counter = 1:length(s(:,1))
    s(counter,1) = s(counter,1) - 4.8032 ;
end

% Final Discarding
s(644:end,:) = [];

for counter1 = 1:length(Vg)
% Bias Condition
OutGenerate.surface_potential =  Vg(counter1);
% Initial Guess Value
[Vbetterguess , A]   = DG_initial(OutGenerate);
% Interpolation 
x_reqd    = 0:0.1:8.7 ;
ab_initio = spline(s(:,1),s(:,counter1 +2) , x_reqd );

% Non - Linear Fitting

  fitting_factor_initial = 5;
  lb=0;
  ub=100;
  options.TolFun=1e-20;
  options.TolX=1e-20;
  options.MaxIter=500;
  options.MaxFunEvals=2000;
  fitting_factor_final = lsqnonlin(@Final_fitting,fitting_factor_initial,lb,ub,options,OutGenerate,ab_initio,Vbetterguess,A);
  fprintf('DG fitting parameter = %f',fitting_factor_final);  


end        



   Schrodinger_Poisson  = open('output/trial/dev1_sp/mat_files/file_0.4868.mat');
% end
  
  
  % Plotting
  S_bias = [0 0.471 0.5026 0.5262 0.546 0.5636 0.5797 0.5946 0.6086 0.62187 0.63453 0.6466];
  Extra_electron_1 = [0, -0.001, -0.002, -0.003, -0.004, -0.005, -0.006, -0.007, -0.008, -0.009, -0.010, -0.011];
  figure(1);
  plot(S_bias,Extra_electron_1);
  hold on;
  
  figure(1);
  plot(Vg,SP_Extra_electron);
  hold on;
  
  fitting_factor_final = linspace(0.1 , 10 , 10 );
  for counter2 = 1:length(fitting_factor_final)
  for counter1 = 1:length(Extra_electron)

    OutGenerate.surface_potential =  Vg(counter1);

    [Vbetterguess, A]   = DG_initial(OutGenerate);

    Density_Gradient = scsdg(fitting_factor_final(counter2),OutGenerate,Vbetterguess,A);

    % Calculaing N chagrge integral using simpsons rule

    x = OutGenerate.x;
    h = (x(end) - x(1))*(0.0529)*(10^-7) / length(x) ;
    Integral = 0;
    for k=1:2:(length(x)-2)
    Integral =  Integral + (h/3) * (Density_Gradient(k) + 4*Density_Gradient(k+1) + Density_Gradient(k+2)) * (1/((5.29*10^-9)^3));
    end
    DG_Extra_electron(counter1) = -Integral*(3.84*3.84*10^-16) ;


  end
  c = {'0.1' '1.2' '2.3' '3.4' '4.5' '5.6' '6.7' '7.8' '8.9' '10' };

  figure(1);
  plot(Vg,DG_Extra_electron);
  xlabel('Gate Voltage(V)');
  ylabel('Extra Electron');
  %legend('abinitio','D.G');
  legend('abinitio','S.P','0.1', '1.2', '2.3', '3.4', '4.5', '5.6', '6.7', '7.8' ,'8.9' ,'10');
  hold on;
  end
  
end