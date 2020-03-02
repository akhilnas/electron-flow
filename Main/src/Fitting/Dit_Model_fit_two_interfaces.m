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


csvmat=csvread('Naomi_Cmat.csv');
Vg=csvmat(30:50,1);
Ctot=csvmat(30:50,2:end); %1MHz CV
f=csvmat(1,2:end);
f(isnan(f))=[];

csvmat=csvread('Naomi_Gwmat.csv');
Gpwtot=csvmat(30:50,2:end); %1MHz CV

csvmat=csvread('Naomi_Gmat.csv');
Gtot=csvmat(30:50,2:end); %1MHz CV

%=============================================

%% Reading Ideal File

reply=input('please provide dielectric thickness (in nm) and oxide capacitance (in F cm^-2)\n','s');
[Tox Cox] =strtok(reply)
Tox=str2num(Tox);
Tox=Tox*18.9;
Cox=str2num(Cox);
Ideal_Dir=uigetdir(pwd);
[OutGenerate HFCV]=Ideal_HFCV(Ideal_Dir, Tox);
[forget trap_loc_index]=min(abs(OutGenerate.x-Tox));
% Vg=HFCV.Vg';36
Vpsi=HFCV.psi_s'/13.6;
% Cox=HFCV.Cox
Cs_tmp=HFCV.Cs';
C=HFCV.C;
psi_s=HFCV.psi_s;
Tempr=OutGenerate.Tempr;
% Tempr=0.0258;

%% Making Ev=0, all energies in eV

Ef_Ev = (OutGenerate.Ef(trap_loc_index)-OutGenerate.Ev(trap_loc_index))*13.6;  
Ec = OutGenerate.Ec(trap_loc_index)*13.6+Ef_Ev;
Ev = zeros(1,length(Ec));
Eg = Ec-Ev;

%% Simulation Parameters

Omega     = 2*3.142*f;

Vt        = 1*Tempr*13.6; % thermal energy
Efermi    = psi_s+Ef_Ev;

E=linspace(-0.4,1,500);


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
Fitting_Details.Dit_N_eq=4;
Fitting_Details.sigma_N_eq=1;
Fitting_Details.N_interfaces=2;

fprintf('Generating Initial Guess\n');
[Dit_param0,ub,lb]=fitting_inital_guess_two_interfaces(Fitting_Details,E,Eg);
% Dit_param0(18)=0.1;
% Dit_param0(18)=-0.64;
% Dit_param0(12)=-0.46;
% Dit_param0(14)=-0.64;




%======================================================================================

options.TolFun=1e-6;
options.TolX=1e-6;
options.MaxFunEvals=8000;
options.MaxIter = 500;
options.Display='iter';

figure(1);
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 6 6]);

Dit_param1=lsqnonlin(@Dit_Model_N_interfaces,Dit_param0,lb,ub,options,simulation_param,Fitting_Details,HFCV,Ctot,Gpwtot);
% 
% 
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

% [Vg1,Ctot1,Gtot1,Gpwtot1,Cp1,Gp1,Gpw1]=Dit_admittance_p(HFCV,Simulation_INP);
% save('Dit_fitted','Vg1','Ctot1','Gtot1','Gpwtot1','Cp1','Gp1','Gpw1','E','Dit_param','Dit','sigma_p','f');
% save('Measured','Vg','Ctot','Gtot','Gpwtot');
% 


%% Admittance Calculations

for interface=1:Fitting_Details.N_interfaces
    start_ind=(interface-1)*length(Dit_param)/Fitting_Details.N_interfaces+1;
    end_ind=interface*length(Dit_param)/Fitting_Details.N_interfaces;
    
    [Dit,sigma_p]=Dit_sigma_gen(Fitting_Details,Dit_param(start_ind:end_ind),E,Vt);
    DDit{interface}=Dit;
    ssigma_p{interface}=sigma_p;

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

Ys{interface}=Dit_admittance_p_N_interfaces(HFCV,Simulation_INP);
end

Omega=2*pi*f;
Ys_tot=zeros(length(f),length(Vg));
for interface=1:Fitting_Details.N_interfaces
    Ys_tot=Ys_tot+Ys{interface};
end


for jj=1:length(f)
    imw=1i*Omega(jj);

    Ys_tot_tmp=imw*HFCV.Cs'+Ys_tot(jj,:);
%     
    Cp(jj,:)=imag(Ys_tot_tmp)/Omega(jj);
    Gp(jj,:)=imag(Ys_tot_tmp)/Omega(jj);
    Gpw(jj,:)=real(Ys_tot_tmp)/Omega(jj);
    
    Ytot=(1./(1/(imw*Cox)+1./(Ys_tot_tmp)));
    Ytot=1./(1./(Ytot));
    %(1/(Gs+imw*Cser)
    
  
    Ctot1(jj,:)=imag(Ytot)/Omega(jj);
    Gpwtot1(jj,:)=real(Ytot)/Omega(jj);
    Gtot1(jj,:)=real(Ytot);
end

% 
