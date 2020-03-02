% Newton Raphson for Cox estimation
% Status - departed
% Please refer CoxEstimation1.m
% Maintained by : Dhirendra (dhirendra22121987@gmail.com)

clc;
clear all;
% close all;
clear functions;
addpath(genpath(pwd));

Q=1.4142; eps0=0.0796; Vt=0.0258/13.6;

WORK_DIR=pwd;
OutGenerate=generate(WORK_DIR);


H=getH(OutGenerate);
H=sparse(H);


ToxOld=1;
% ToxNew=1.1;
% [Forget interace_index]=min(abs(OutGenerate.x-ToxOld*18.9));
% NTox=interace_index-1;
% xTox=ToxNew/NTox;
% OutGenerate.x=OutGenerate.x+(ToxNew-ToxOld)*18.9;
% OutGenerate.x(1:NTox+1)=linspace(0,ToxNew*18.9,NTox+1);
Cox=4e-6;
check=1;
Maxit=50;
hold on;
CoxSim=getCox(H,OutGenerate)
iter=0;
while(check)
    iter=iter+1;
ToxNew=ToxOld*1e-7+(CoxSim-Cox)*((ToxOld*1e-7)^2)/(3.9*8.85e-14);
ToxNew=ToxNew*1e7;
% ToxNew=round(ToxNew*1e4)/1e4;
[Forget interace_index]=min(abs(OutGenerate.x-ToxOld*18.9))
NTox=interace_index;
OutGenerate.x=OutGenerate.x+(ToxNew-ToxOld)*18.9;
OutGenerate.x(1:NTox)=linspace(0,ToxNew*18.9,NTox);
H=getH(OutGenerate);
H=sparse(H);
clear NRscs3;
clear functions;
CoxSim=getCox(H,OutGenerate);

fprintf('Iteration no %d    Tox=%f   Cox=%f ================',iter,ToxNew,CoxSim*1e6);
input('')

ToxOld=ToxNew;

if iter>500
    check=0;
end

if abs(CoxSim-Cox)<1e-10
    check=0;
end

end


% while(check)
%     
%     iter=iter+1;
%     CoxSim=getCox(H,OutGenerate) % Cox from simulator
%     Tox
%     input('');
%     Tox_tmp=Tox;
%     Tox=Tox*1e-7+(CoxSim-Cox)*((Tox*1e-7)^2)/(3.9*8.85e-14);
%     Tox=Tox*1e7;
%     [Forget interace_index]=min(abs(OutGenerate.x-Tox_tmp));
%     NTox=interace_index;
%     OutGenerate.x=OutGenerate.x+(Tox-Tox_tmp)*18.9;
%     OutGenerate.x(1:NTox)=linspace(0,Tox*18.9,NTox);
%     H=getH(OutGenerate);
%     H=sparse(H);
%     
% %     iter=iter+1;
% %     Tox
% %     CoxSim=getCox(H,OutGenerate) % Cox from simulator
% % %     figure(1);
% %     plot(iter,Tox,'o');
% %     Tox_tmp=Tox;
% %     Tox=Tox*1e-7+(CoxSim-Cox)*((Tox*1e-7)^2)/(3.9*8.85e-14);
% % %     Tox=Tox*1e-7+(CoxSim-Cox)/dCTox;
% %     Tox=Tox*1e7;
% % %     input('','s');
% %     [Forget interace_index]=min(abs(OutGenerate.x-Tox_tmp));
% %     NTox=interace_index-1;
% %     xTox=Tox/NTox;
% %     OutGenerate.x=OutGenerate.x+(Tox-Tox_tmp)*18.9;
% %     OutGenerate.x(1:NTox+1)=linspace(0,Tox*18.9,NTox+1);
% %     H=getH(OutGenerate);
% %     H=sparse(H);
%     
%     fprintf('Iteration no ============================= %d\n',iter);
%     
%     if Maxit>500
%         check=0;
%     end
%     if abs(CoxSim-Cox)<1e-8
%         check=0;
%     end
%     
% end

% 
% 
% if ( OutGenerate.control.Solvers.SP || OutGenerate.control.Solvers.P )
%     for i = 1:size(OutGenerate.surface.bias,2)
%         try
%             
%             NRscs3(H,OutGenerate,OutGenerate.surface.bias{i});
%         catch err
%             idSegLast = regexp(err.identifier, '(?<=:)\w+$', 'match');
%             if  strcmp(idSegLast,'ConvergenceFailed')
%                 fprintf('\n'); fprintf(err.message); fprintf('\n\n');
%                 return
%             else
%                 rethrow(err);
%             end
%         end
%     end
% 
% output_dir  = OutGenerate.control.Output_directory; 
% file_prefix = OutGenerate.control.file_prefix;
% 
% fprintf('\nSaving the results in:\n %s \n',output_dir);
% Generate_text_CV(output_dir,file_prefix);
% generate_text_output(output_dir);
% 
% end
% 
% 
% if OutGenerate.control.Solvers.Trap
%     
%     output_dir = OutGenerate.control.Output_directory;
%     files  = dir(fullfile(output_dir, 'mat_files','*.mat'));
%     
%     for i=1:length(files)
%         input_file = fullfile(output_dir, 'mat_files',files(i).name);
%         load(input_file);        
%     end
% end
% 
% fprintf('\n\nSimulation completed...\n');