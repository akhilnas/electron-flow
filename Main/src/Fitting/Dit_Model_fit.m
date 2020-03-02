clc;
clear all;
close all;

%% Global constants
Kb = 1; Q_SI = 1.6e-19; Ab_CGS = 0.529e-8; Q=sqrt(2);

%% Reading constraints

if ~exist('pathname')
    pathname=pwd;
else
    if pathname==0
        pathname=pwd;
    end
end



Vg_start=-0.1;
Vg_end=0.85;

csvmat=csvread('CV_Data/synth_S17/Cmat.csv');
Vg_tmp=csvmat(:,1);
[~, Vg_start_Ind]=min(abs(Vg_tmp-Vg_start));
[~, Vg_end_Ind]=min(abs(Vg_tmp-Vg_end));

% csvmat=csvread('CV_Data/Naomi_Cmat_synth_S17.csv');
Vg=csvmat(Vg_start_Ind:Vg_end_Ind,1);
Ctot=csvmat(Vg_start_Ind:Vg_end_Ind,2:end-1); %1MHz CV
f=csvmat(1,2:end-1);
f(isnan(f))=[];

% csvmat1=csvread('synth_SISPAD_stretchout.csv');
% C1M=csvmat1(Vg_start_Ind-1:Vg_end_Ind-1,end); %1MHz CV
C1M=csvmat(Vg_start_Ind:Vg_end_Ind,end-1); %1MHz CV
C10M=csvmat(Vg_start_Ind:Vg_end_Ind,end); %1MHz CV

csvmat2=csvread('CV_Data/synth_S17/Gwmat.csv');
Gpwtot=csvmat2(Vg_start_Ind:Vg_end_Ind,2:end-1); %1MHz CV

csvmat3=csvread('CV_Data/synth_S17/Gmat.csv');
Gtot=csvmat3(Vg_start_Ind:Vg_end_Ind,2:end-1); %1MHz CV

%% Reading Ideal File

reply=input('please provide dielectric thickness (in nm) and oxide capacitance (in F cm^-2)\n','s');
[Tox Cox] =strtok(reply)
Tox=str2num(Tox);
Tox=Tox*18.9;
Cox=str2num(Cox);


% Ideal_Dir=uigetdir(pwd);
if (exist('Ideal_Dir')==0)
    Ideal_Dir=uigetdir(pwd);
end
if (exist('Ideal_Dir')==1)
    if (Ideal_Dir==0)
        Ideal_Dir=uigetdir(pwd);
    end
end

[OutGenerate HFCV_main]=Ideal_HFCV(Ideal_Dir, Tox);
[forget trap_loc_index]=min(abs(OutGenerate.x-Tox));
% Vg=HFCV.Vg';36
Vpsi=HFCV_main.psi_s'/13.6;
% Cox=HFCV.Cox
Cs_tmp=HFCV_main.Cs';
C=HFCV_main.C;
psi_s=HFCV_main.psi_s;
Tempr=OutGenerate.Tempr;
% Tempr=0.0258;


Vg_interpolated=interp1(HFCV_main.C,HFCV_main.Vg,C1M);

HFCV=get_HFCV_interpolated(Vg_interpolated,HFCV_main);
 
%% Making Ev=0, all energies in eV

Ef_Ev = (OutGenerate.Ef(trap_loc_index)-OutGenerate.Ev(trap_loc_index))*13.6;  
Ec = OutGenerate.Ec(trap_loc_index)*13.6+Ef_Ev;
Ev = zeros(1,length(Ec));
Eg = Ec-Ev;

%% Simulation Parameters

Omega     = 2*3.142*f;

Vt        = 1*Tempr*13.6; % thermal energy
Efermi    = psi_s+Ef_Ev;

E=linspace(-0.4,1.12,500);


%% Model fitting

AAA=[Vg Ctot];
AAA=sortrows(AAA,1);

Ctot=AAA(:,2:end)';

AAA=[Vg Gpwtot];
AAA=sortrows(AAA,1);

Gpwtot=AAA(:,2:end)';


AAA=[Vg Gtot];
AAA=sortrows(AAA,1);

Gtot=AAA(:,2:end)';


clr=[rand(length(Vg),1) rand(length(Vg),1) rand(length(Vg),1)];

simulation_param.f=f;
simulation_param.Omega=Omega;
simulation_param.E=E;
simulation_param.Efermi=Efermi;
simulation_param.Vt=Vt;
simulation_param.Vg=HFCV.Vg;
simulation_param.HFCV=HFCV;
simulation_param.Cox=Cox;
simulation_param.Vm=Vg;
simulation_param.clr=clr;


Fitting_Details.Dit_fit_eq='Gaussian';
Fitting_Details.sigma_fit_eq='Exponential';
Fitting_Details.Dit_N_eq=2;
Fitting_Details.sigma_N_eq=1;

fprintf('Generating Initial Guess\n');
[Dit_param0,ub,lb]=fitting_inital_guess(Fitting_Details,E,Eg);
% Dit_param0(9)=0;
% Dit_param0(8)=-16;
% Dit_param0(15)=0;
% ub(15)=0.0000001;
% lb(15)=-0.0000001;
% Dit_param0(12)=-0.46;
% Dit_param0(14)=-0.64;




%======================================================================================

options.TolFun=1e-12;
options.TolX=1e-12;
options.MaxFunEvals=12000;
options.MaxIter = 500;
options.Display='iter';

figure(1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 6 6]);

Dit_param1=lsqnonlin(@Dit_Model,Dit_param0,lb,ub,options,simulation_param,Fitting_Details,HFCV,Ctot,Gpwtot,Gtot);


Dit_param=Dit_param1;


%% Dit and Capture Cross Section

[Dit,sigma_p]=Dit_sigma_gen(Fitting_Details,Dit_param,E,Vt);
sigma_n=1e-16*ones(1,length(E));
vthn=2e7*ones(1,length(E));
vthp=2e7*ones(1,length(E));

tau_n0=1./(sigma_n.*vthn);
tau_p0=1./(sigma_p.*vthp);

%% Admittance calculations
Simulation_INP.f=f;
Simulation_INP.E=E;
Simulation_INP.Dit=Dit;
Simulation_INP.tau_n0=tau_n0;
Simulation_INP.tau_p0=tau_p0;
Simulation_INP.Cox=Cox;

[Vg1,Ctot1,Gtot1,Gpwtot1,Cp1,Gp1,Gpw1,Cit,Git,Gpwit]=Dit_admittance_p1(HFCV,Simulation_INP);


save('Dit_profiles/Dit_synth_S17_GVfit.mat','Vg1','Ctot1','Gtot1','Gpwtot1','Cp1','Gp1','Gpw1','E','Dit_param','Dit','sigma_p','f');
save('Measured_synth_S17','Vg','Ctot','Gtot','Gpwtot');


csvmat=csvread('CV_Data/synth_S17/Cmat1.csv');
Vg_tmp=csvmat(:,1);
[~, Vg_start_Ind]=min(abs(Vg_tmp-Vg_start));
[~, Vg_end_Ind]=min(abs(Vg_tmp-Vg_end));
% Vg=csvmat(Vg_start_Ind:Vg_end_Ind,1);
Csemi=csvmat(Vg_start_Ind:Vg_end_Ind,2:end-1); 
% C1Msemi=csvmat(Vg_start_Ind:Vg_end_Ind,end-1); %1MHz CV
% 
% Cstr_out_tmp=C1Msemi-Cit(:,end);
% Csemi=csvmat(Vg_start_Ind:Vg_end_Ind,2:end-1);
Csideal=Csemi-Cit';
Cs_str_out=Csideal(:,end);
C_str_out=Cox*Cs_str_out./(Cs_str_out+Cox);



[OutGenerate HFCV_main]=Ideal_HFCV(Ideal_Dir, Tox);
[forget trap_loc_index]=min(abs(OutGenerate.x-Tox));
% Vg=HFCV.Vg';36
Vpsi=HFCV_main.psi_s'/13.6;
% Cox=HFCV.Cox
Cs_tmp=HFCV_main.Cs';
C=HFCV_main.C;
psi_s=HFCV_main.psi_s;
Tempr=OutGenerate.Tempr;
% Tempr=0.0258;


Vg_interpolated=interp1(HFCV_main.C,HFCV_main.Vg,C_str_out);

HFCV=get_HFCV_interpolated(Vg_interpolated,HFCV_main);
 
%% Making Ev=0, all energies in eV

Ef_Ev = (OutGenerate.Ef(trap_loc_index)-OutGenerate.Ev(trap_loc_index))*13.6;  
Ec = OutGenerate.Ec(trap_loc_index)*13.6+Ef_Ev;
Ev = zeros(1,length(Ec));
Eg = Ec-Ev;

%% Simulation Parameters

Omega     = 2*3.142*f;

Vt        = 1*Tempr*13.6; % thermal energy
Efermi    = psi_s+Ef_Ev;

E=linspace(-0.4,1.12,500);


%% Model fitting

% AAA=[Vg Ctot];
% AAA=sortrows(AAA,1);
% 
% Ctot=AAA(:,2:end)';
% 
% AAA=[Vg Gpwtot];
% AAA=sortrows(AAA,1);
% 
% Gpwtot=AAA(:,2:end)';
% 
% 
% AAA=[Vg Gtot];
% AAA=sortrows(AAA,1);
% 
% Gtot=AAA(:,2:end)';


clr=[rand(length(Vg),1) rand(length(Vg),1) rand(length(Vg),1)];

simulation_param.f=f;
simulation_param.Omega=Omega;
simulation_param.E=E;
simulation_param.Efermi=Efermi;
simulation_param.Vt=Vt;
simulation_param.Vg=HFCV.Vg;
simulation_param.HFCV=HFCV;
simulation_param.Cox=Cox;
simulation_param.Vm=Vg;
simulation_param.clr=clr;


Fitting_Details.Dit_fit_eq='Gaussian';
Fitting_Details.sigma_fit_eq='Exponential';
Fitting_Details.Dit_N_eq=2;
Fitting_Details.sigma_N_eq=1;

fprintf('Generating Initial Guess\n');
[Dit_param0,ub,lb]=fitting_inital_guess(Fitting_Details,E,Eg);
% Dit_param0(9)=0;
% Dit_param0(8)=-16;
% Dit_param0(15)=0;
% ub(15)=0.0000001;
% lb(15)=-0.0000001;
% Dit_param0(12)=-0.46;
% Dit_param0(14)=-0.64;




%======================================================================================

options.TolFun=1e-12;
options.TolX=1e-12;
options.MaxFunEvals=12000;
options.MaxIter = 500;
options.Display='iter';

figure(1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 6 6]);

Dit_param1=lsqnonlin(@Dit_Model1,Dit_param0,lb,ub,options,simulation_param,Fitting_Details,HFCV,Ctot,Gpwtot,Gtot);


Dit_param=Dit_param1;


%% Dit and Capture Cross Section

[Dit,sigma_p]=Dit_sigma_gen(Fitting_Details,Dit_param,E,Vt);
sigma_n=1e-16*ones(1,length(E));
vthn=2e7*ones(1,length(E));
vthp=2e7*ones(1,length(E));

tau_n0=1./(sigma_n.*vthn);
tau_p0=1./(sigma_p.*vthp);

%% Admittance calculations
Simulation_INP.f=f;
Simulation_INP.E=E;
Simulation_INP.Dit=Dit;
Simulation_INP.tau_n0=tau_n0;
Simulation_INP.tau_p0=tau_p0;
Simulation_INP.Cox=Cox;

[Vg2,Ctot2,Gtot2,Gpwtot2,Cp2,Gp2,Gpw2,Cit2,Git2,Gpwit2]=Dit_admittance_p1(HFCV,Simulation_INP);