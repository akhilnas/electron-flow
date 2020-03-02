%% Single level interface trap Admittance
%
%
%
%

%%

clc;
clear all;
% close all;

%% Global constants
Kb = 1; Q_SI = 1.60217e-19; Ab_CGS = 0.529e-8;

%% Dit miniparser

if (exist('Ideal_Dir')==0)
    Ideal_Dir=uigetdir(pwd);
end
if (exist('Ideal_Dir')==1)
    if (Ideal_Dir==0)
        Ideal_Dir=uigetdir(pwd);
    end
end
reply=input('please provide Tox (in nm) and Cox (in F cm^-2)\n','s');
[Tox Cox] =strtok(reply)
Tox=str2num(Tox);
Tox=Tox*18.9;
Cox=str2num(Cox);

[OutGenerate HFCV]=Ideal_HFCV(Ideal_Dir, Tox);
[forget trap_loc_index]=min(abs(OutGenerate.x-Tox));
Vg=HFCV.Vg';
Vpsi=HFCV.psi_s'/13.6;
Cox=HFCV.Cox;
Cs_tmp=HFCV.Cs';
C=HFCV.C;
Tempr=OutGenerate.Tempr;

%% Simulation parameters
V_applied = 0.4247;              % fermi level bias
f         = logspace(-1,10,100);
Omega     = 2*3.1425*f;
Cs=interp1(Vg,Cs_tmp,V_applied);

%% Interface trap parameters

Tempr  = (0.0019/300)*100;
E_Ev   = linspace(0,0.67,100);                              % energy spectrum
Ec     = (OutGenerate.Ec(trap_loc_index)-OutGenerate.Ev(trap_loc_index))*13.6;              % Ec w.r.t Ev
sigma  = 1e-16;                                             % capture cross section
vth    = 2e7;                                               % Thermal velocity cm/s
Neff   = OutGenerate.Nv(end)*((0.0529e-7)^(-3));
tau_p  = (vth*sigma*Neff*exp(-E_Ev/(Kb*Tempr*13.6))).^(-1);
tau_n  = (vth*sigma*Neff*exp(-(Ec-E_Ev)/(Kb*Tempr*13.6))).^(-1);
tau_E(1:length(E_Ev/2))=tau_p(1:length(E_Ev/2));
tau_E(length(E_Ev)/2+1:length(E_Ev))=tau_n(length(E_Ev)/2+1:length(E_Ev));



E      = 0.1;                                              % w.r.t Ev  (E-Ev)
Dit    = 1e13;                                              % in eV^-1 cm^-2
tau_it = interp1(E_Ev,tau_E,E);

%% C and G/w calculations

imw    = 1i*Omega;
Cit    = Q_SI*Dit./(1+(Omega*tau_it).^2);
Git_w  = Q_SI*Dit*(Omega*tau_it)./(1+(Omega*tau_it).^2);
% Yit    = imw.*Cit + Git_w .* Omega;
% Ys     = imw.*Cs+Yit;
% Ytot   = 1/Cox+1./Ys;
% 
% Ctot   = real(Ytot);
% Gtot   = imag(Ytot);
% Gtot_w = Gtot./Omega;

%% plotting
figure(1);
set(0,'DefaultAxesFontSize',16);

subplot(1,2,1);
semilogx(f,Git_w);
xlabel('frequency (Hz)'); hold all;
ylabel('G_{it}/\omega (F cm^{-2})');

subplot(1,2,2);
semilogy(E_Ev,tau_E);
xlabel('E-Ev (eV)');
ylabel('\tau_{it} (sec)');
xlim([E_Ev(1) E_Ev(end)]);





