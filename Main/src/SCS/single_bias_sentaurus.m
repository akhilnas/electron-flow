% Newton - Raphson implementation of CVsimulator
% Started on 21 May 2014 by dhirendra

% Maintained by : Dhirendra (dhirendra22121987@gmail.com)

function single_bias_abinitio(input_file)

OutGenerate=generate(input_file);


output_dir  = OutGenerate.control.Output_directory; 
file_prefix = OutGenerate.control.file_prefix;



  [Vbetterguess , A]   = DG_initial(OutGenerate);
  %ab_initio  = open('output/trial/dev1_sp/mat_files/file_0.1.mat');
  %fid = fopen('Akhil/Si/eDensity_5nm.csv');
  %mydata = textscan(fid,'%s %s %s %s %s %s %s %s %s %s','Delimiter',',','CollectOutput',1);
  %mydata{1,1}(1,:) = [];
  inputdata = csvread('Akhil/Micrometre/Si/eD_10nm.csv',1);
  inputdata2 = csvread('Akhil/Micrometre/Si/EP_10nm.csv',1);
  ndist = inputdata(:,3);
  xdist = inputdata(:,1);
  vdist = inputdata2(:,3);
  xdistv = inputdata2(:,1);
  format long,xdist,inputdata,vdist;
  % Recorrecting for channel width
  x_input = xdist(14:38);
  n_charge_input = ndist(14:38);
  vinput = vdist(14:38);
  xinputv = xdistv(14:38);
  % Correcting
  x_input(2) = [];
  xinputv(2) = [];
  n_charge_input(2) = [];
  vinput(2) = [];
  x_input(end-1) = [];
  xinputv(end-1) = [];
  n_charge_input(end-1) = [];
  vinput(end-1) = [];
  % Recorrecting x coordinates
  Factor = -x_input(1)*ones(length(x_input),1);
  x_input = x_input+ Factor;
  x_input = x_input.*1000;
  x_input = x_input./0.0529;
  n_charge_input = n_charge_input.*((5.29*10^-9)^3);
  % Recorrecting x coordinates
  Factor1 = -xinputv(1)*ones(length(xinputv),1);
  xinputv = xinputv+ Factor1;
  xinputv = xinputv.*1000;
  figure(1);
  plot(x_input,n_charge_input/((5.29*10^-9)^3));
  figure(2);
  plot(xinputv,vinput);
  % Performing Interpolation/Extrapolation
  % Condition Check if Interpolation or Extrapolation is to be used
       n_charge_ref = interp1(x_input,n_charge_input,OutGenerate.x,'linear');
  % Correcting
        n_charge_ref(end) =0;
  fitting_factor_initial = 10;
  lb=-100;
  ub=300;
  options.TolFun=1e-40;
  options.TolX=1e-50;
  options.MaxIter=100;
  options.MaxFunEvals=1000;
  fitting_factor_final = lsqnonlin(@DG_fitting_abinitio,fitting_factor_initial,lb,ub,options,OutGenerate,n_charge_ref,Vbetterguess,A);
  fprintf('DG fitting parameter = %f',fitting_factor_final);
figure(1);
plot(OutGenerate.x*0.0529,n_charge_ref/((5.29*10^-9)^3),'--','LineWidth',2);
ylabel('Electron Density(cm-3)');
xlabel('x(nm)');
hold all;
scsdg(fitting_factor_final,OutGenerate,Vbetterguess,A);

fprintf('\n\nSimulation completed...\n');
DG = open('output/trial/dev1_sp/mat_files/file_0_DG.mat');
figure(1);
plot(DG.x*0.0529,DG.n_charge_DG/((5.29*10^-9)^3),'LineWidth',2);
legend('S.P','D.G');
hold all;
 [ VDG ] = Poisson(DG.n_charge_DG,DG.p_charge_DG,OutGenerate);
 figure(3);
 plot(xinputv,vinput,'--','LineWidth',2);
 ylabel('Potential(Volts)');
 xlabel('x(nm)');
 %axis([0 5 -0.2 0]);
 hold all;
 figure(3);
 plot(DG.x*0.0529,VDG*19.2,'LineWidth',2);
 legend('S.P','D.G');
 hold all;
end
