function diff1=Dit_Model(Dit_param,simulation_param,Fitting_Details,HFCV,Cmeas,Gpwmeas,Gmeas)

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

[Vg,Ctot,Gtot,Gpwtot,Cp,Gp,Gpw]=Dit_admittance_p(HFCV,Simulation_INP);



% size(Ctot)
% size(Cmeas)

% dd=diff(Vg)./diff(Vm);
% dd(end+1)=dd(end);
diff1=1e6*(Ctot-Cmeas);

% diff1=[diff1; 1e2*(dd'-1)];

abcd=0;
if (abcd)
figure(1)
subplot(3,2,1);
hold off;
for i=1:length(Vg)
    semilogx(f,Ctot(:,i),'-s','color',clr(i,:)); hold all;
    semilogx(f,Cmeas(:,i),'-o','color',clr(i,:)); hold all;
end
xlabel('frequency (Hz)');
ylabel('Capacitance (F/cm^2)')


% figure(2)
subplot(3,2,2);
hold off;
for i=1:length(Vg)
    semilogx(f,Gpwtot(:,i),'-s','color',clr(i,:)); hold all;
    semilogx(f,Gpwmeas(:,i),'-o','color',clr(i,:)); hold all;
end
xlabel('frequency (Hz)');
ylabel('Conductance (S/cm^2)')

% figure(3)
subplot(3,2,3);
semilogy(E,Dit);
xlabel('E-Ev (eV)');
ylabel('Dit /eV/cm^2');

% figure(4) 
subplot(3,2,4);
semilogy(E,sigma_p);
xlabel('E-Ev (eV)');
ylabel('Capture Cross Section (cm^2)');
% figure(5)
subplot(3,2,5);
plot(Vg,Ctot,'-',Vm,Cmeas,'o');
xlabel('gate bias (V)');
ylabel('Capacitance (F/cm^2)');
subplot(3,2,6);
plot(Vg,Gtot,'-',Vm,Gmeas,'o');
xlabel('gate bias (V)');
ylabel('Conductance (S/cm^2)');
end

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