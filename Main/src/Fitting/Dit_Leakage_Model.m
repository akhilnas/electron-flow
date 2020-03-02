function diff=Model4(Dit_param,simulation_param,HFCV,Cmeas,Gpwmeas)

%% Global constants
Kb = 1; Q_SI = 1.6e-19; Ab_CGS = 0.529e-8; Q=sqrt(2);

%% Simulation Parameters
f=simulation_param.f;
Omega=simulation_param.Omega;
E=simulation_param.E;
Efermi=simulation_param.Efermi;
Vg1=simulation_param.Vg;
Vt=simulation_param.Vt;
HFCV=simulation_param.HFCV;
Cox=simulation_param.Cox;
clr=simulation_param.clr;
Vm=simulation_param.Vm;

%% Dit and Capture Cross Section

sigma_n=1e-16*ones(1,length(E));
vthn=2e7*ones(1,length(E));
vthp=2e7*ones(1,length(E));


Dit=10^Dit_param(1)*exp(-(E-Dit_param(2)).^2/(Dit_param(3)))+10^Dit_param(4)*exp(-(E-Dit_param(5)).^2/(Dit_param(6)))+10^Dit_param(7)*exp(-(E-Dit_param(8)).^2/(Dit_param(9)));
Dit=Dit+10^Dit_param(10);

% Dit=Dit_param(1)*exp(-(E-Dit_param(2)).^2/(Dit_param(3)))+Dit_param(4)*exp(-(E-Dit_param(5)).^2/(Dit_param(6)))+Dit_param(7)*exp(-(E-Dit_param(8)).^2/(Dit_param(9)));
% Dit=Dit+Dit_param(10);

% % sigma_p=2e-18*exp(-E*(-0.46)/Vt)+1e-18;
% sigma_p=10^Dit_param(11)*exp(-E*(Dit_param(12))/Vt)+10^Dit_param(13);
% sigma_p=10^Dit_param(11)*exp(-E*(Dit_param(12))/Vt)+10^Dit_param(13);
% sigma_p=10^Dit_param(11)*ones(1,length(E));
% 
%=======================================================
% %gaussian fitting sigma
sigma_p=10^Dit_param(11)*exp(-(E-Dit_param(12)).^2/(Dit_param(13)))+10^Dit_param(14)*exp(-(E-Dit_param(15)).^2/(Dit_param(16)))+10^Dit_param(17)*exp(-(E-Dit_param(18)).^2/(Dit_param(19)));
sigma_p=sigma_p+10^Dit_param(20);
%==========================================================

tau_n0=1./(sigma_n.*vthn);
tau_p0=1./(sigma_p.*vthp);

%% Admittance calculations
Simulation_INP.f=f;
Simulation_INP.E=E;
Simulation_INP.Dit=Dit;
Simulation_INP.tau_n0=tau_n0;
Simulation_INP.tau_p0=tau_p0;
Simulation_INP.Cox=Cox;
% Simulation_INP.sigma_G=Dit_param(21:end);
% Simulation_INP.Gs=Dit_param(21);
% Simulation_INP.Cser=Dit_param(22);

[Vg,Ctot,Gtot,Gpwtot,Cp,Gp,Gpw]=Dit_admittance_p_leakage(HFCV,Simulation_INP);

% size(Ctot)
% size(Cmeas)

diff=1e6*(Gpwtot-Gpwmeas)+1e6*(Ctot-Cmeas);

figure(1)
subplot(3,2,1);
hold off;
for i=1:length(Vg)
    semilogx(f,Ctot(:,i),'-s','color',clr(i,:)); hold all;
    semilogx(f,Cmeas(:,i),'-o','color',clr(i,:)); hold all;
end


% figure(2)
subplot(3,2,2);
hold off;
for i=1:length(Vg)
    semilogx(f,Gpwtot(:,i),'-s','color',clr(i,:)); hold all;
    semilogx(f,Gpwmeas(:,i),'-o','color',clr(i,:)); hold all;
end

% figure(3)
subplot(3,2,3);
semilogy(E,Dit);

% figure(4)
subplot(3,2,4);
semilogy(E,sigma_p);

% figure(5)
subplot(3,2,5);
plot(Vg,Ctot,'-',Vm,Cmeas,'o');

% subplot(3,2,6);
% plot(Vg,Gmeas-Gtot);
% 
% figure(2)
% semilogx(f,Gpwmeas-Gpwtot);
% 
% figure(3)
% semilogx(f,Gmeas-Gtot)
% 
% Gdiff=Gmeas-Gtot;
% 
% figure(1);semilogx(f,Ctot,'-s','color',hsv(23),f,Cmeas,'-o','color',hsv(23))
% figure(2);semilogx(f,Gpwtot,'-s','color',hsv(23),f,Gpwmeas,'-o','color',hsv(23))
% % figure(3);plot(Vg,Ctot,'-',Vg,Cmeas,'o');
% figure(3);semilogy(E,sigma_p);
% figure(4);semilogy(E,Dit);

% figure(5);plot(HFCV.Vg,Ctot)
end