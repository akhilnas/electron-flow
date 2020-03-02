% Newton's Method for Cox estimation
% Status - Ready
% In some cases Newton Method may diverge
% Currently on a single dielectric thickness is implemented. Need to work for inter layer thickness estimation 
%
% Created by : Dhirendra (dhirendra22121987@gmail.com)
% Maintained by : Dhirendra (dhirendra22121987@gmail.com)

clc;
clear all;
% close all;
clear functions;
addpath(genpath(pwd));

%% Mathematical and Physical constants
Q=1.4142; eps0=0.0796; Vt=0.0258/13.6;

%% Accepting Inputs
fprintf('Please provide following inputs as a space separated list\n');
fprintf('Make sure Tox in intial guess and Tox in initial input file are same\n');
fprintf('Initial Guess for Tox (in nm), Cox in F/(cm^-2), Cox error Tolerance in F/(cm^-2), Maximum no of Iterations\n');
fprintf(' Example : 1 3e-6 1e-10 50\n');
reply=input('','s');
[Tox remain]=strtok(reply);
[Cox remain]=strtok(remain);
[Tolerance Maxit]=strtok(remain);
ToxOld=eval(Tox);
Cox=eval(Cox);
Tolerance=eval(Tolerance);
Maxit=eval(Maxit);
WORK_DIR=pwd;
OutGenerate=generate(WORK_DIR);


%% Start
H=getH(OutGenerate);
H=sparse(H);

check=1;
CoxSim=getCox(H,OutGenerate);
iter=0;

while(check)
    iter=iter+1;
    ToxNew=ToxOld*1e-7+(CoxSim-Cox)*((ToxOld*1e-7)^2)/(3.9*8.85e-14); % Newton Raphson Method C=eps/tox. dC/dtox=-eps/(tox^2)
    ToxNew=ToxNew*1e7; % Tox in nm
    % ToxNew=round(ToxNew*1e4)/1e4;
    [Forget interace_index]=min(abs(OutGenerate.x-ToxOld*18.9)); % finding interface location and so the no of points in the dielectric
    NTox=interace_index;
    OutGenerate.x=OutGenerate.x+(ToxNew-ToxOld)*18.9; % Modifying the grid
    OutGenerate.x(1:NTox)=linspace(0,ToxNew*18.9,NTox);
    H=getH(OutGenerate);
    H=sparse(H);
    CoxSim=getCox(H,OutGenerate);
    
    fprintf('Iteration no %d    Tox=%f   Cox=%f ================\n',iter,ToxNew,CoxSim*1e6);
    
    ToxOld=ToxNew;
    
    if iter>Maxit
        check=0;
    end
    
    if abs(CoxSim-Cox)<Tolerance
        check=0;
    end
    
end