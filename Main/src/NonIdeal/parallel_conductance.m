% This file is written to examine the Dit Admittance of MOS.
% As different from previous script in simulator_v2, this script attempts
% to make the code parser free. Inputs are accepted interactively.
%
%
% Created by : Dhirendra (dhirendra22121987@gmail.com)
% Maintained by : Dhirendra (dhirendra22121987@gmail.com)

clc;
% clear all;
close all;

%% Global constants
Kb = 1; Q_SI = 1.60217e-19; Ab_CGS = 0.529e-8;

%% Dit miniparser

% fprintf('Please provide, \n IdealCV directory, \n dielectric thickness in nm, \n Cox in F cm^{-2}, \n\n example: \n /home/dhirendra/simulator/dev_lab/device1 1 3.4515e-6\n \n');
% reply=input('','s');

% [Ideal_Dir reply] =strtok(reply);
if (exist('Ideal_Dir')==0)
    Ideal_Dir=uigetdir(pwd);
end
if (exist('Ideal_Dir')==1)
    if (Ideal_Dir==0)
        Ideal_Dir=uigetdir(pwd);
    end
end
reply=input('please provide dielectric thickness (in nm) and oxide capacitance (in F cm^-2)\n','s');
[Tox Cox] =strtok(reply)
Tox=str2num(Tox);
Tox=Tox*18.9;
Cox=str2num(Cox);

[OutGenerate HFCV]=Ideal_HFCV(Ideal_Dir, Tox);
[forget trap_loc_index]=min(abs(OutGenerate.x-Tox));
Vg=HFCV.Vg';
Vpsi=HFCV.psi_s'/13.6;
% Cox=HFCV.Cox;
Cs_tmp=HFCV.Cs';
C=HFCV.C;
Tempr=OutGenerate.Tempr;


% CVm=importdata('../dev_lab/CV_S5.dat'); % import experimental data
%% Actual script


    % stretchout
    
%     reply =input('Please provide Energy range (in eV w.r.t Ev, Ev=0) \n','s')
    [E_start reply]=strtok(reply);
    [E_end E_step]= strtok(reply);
    E_start=-0;
    E_end=0.67;
    E_step=100;
    psi=linspace(E_start,E_end,E_step);
    zeta=psi/13.6-(OutGenerate.Ei(trap_loc_index)-OutGenerate.Ev(trap_loc_index));
    zeta1=zeta;
    
%     reply = input('Enter Dit profile as a matlab function \n','s');
%     Ditd=1e13*ones(1,length(zeta));
%     Ditd=3e10+0.7e11*normpdf(zeta*13.6,0,0.15);
    
    Ditd = 2.5e11*ones(1,length(zeta));
    Dita=Ditd;
    
    Ei_s = OutGenerate.Ei(trap_loc_index) - Vpsi;
    Ef_s = OutGenerate.Ef(trap_loc_index);
    
    Qit = zeros(1,length(Vpsi));
    
    for i=1:length(Vpsi)
        temp = 1 + exp( (zeta+Ei_s(i)-Ef_s)/(Kb*Tempr) );
        
        Qit(i) = Q_SI * trapz(zeta*13.6,Ditd)  ...
            - Q_SI * trapz(zeta*13.6,(Dita+Ditd)./temp);
    end
    
    Vg_new = Vg - Qit/Cox ;
    Phi_s = Ef_s - Ei_s;
    psi_s=HFCV.psi_s;
    
%     Vm=CVm(:,1);
%     Cm=CVm(:,5);
    
    ZBB=interp1(psi_s,Vg,0);
    if (isnan(ZBB))
        ZBB=-1.5;
    end
    CFB=interp1(Vg,C,ZBB);
%     VFB=interp1(Cm(1:end-10),Vm(1:end-10),CFB*1e6);
    VFB1=interp1(C(1:end),Vg_new(1:end),CFB);
%     Vm=Vm-VFB+ZBB;
    
    
%     figure(1)
%     set(gcf, 'Position', [0 0 550 550])
%     plot(Vg,C,'b','linewidth',2);
%     hold on;
%     plot(Vg_new-VFB1+ZBB,C,'r','linewidth',2)
%     hold on;
% %     plot(Vm,Cm*1e-6,'k');
%     plot(Vm,Cm*1e-6,'sk');
%     hold off;
%     % set(AX1,'fontsize',15);
%     set(gca,'fontsize',20)
%     xlabel('gate bias (V)');
%     ylabel('capacitance (F cm^{-2})');
%     xlim([-2.5 2]);
%     grid on;
%     set(gca,'gridlinestyle','--');
%     legend('Ideal HFCV','Dit stretchout','500 KHz');
    
   


% ==================== frequency dependant admittance =====================

% ==================== Use actual Tau =====================================
Ei_s_Ev_s=OutGenerate.Ei(trap_loc_index)-OutGenerate.Ev(trap_loc_index);
Ec_s_Ei_s=OutGenerate.Ec(trap_loc_index)-OutGenerate.Ei(trap_loc_index);

sigma=1e-12;
Vth=2.3e7;
Nc=OutGenerate.Nc(trap_loc_index)*((0.0529e-7)^(-3));
Nv=OutGenerate.Nv(trap_loc_index)*((0.0529e-7)^(-3));
% tau_p1=(1/(sigma*Vth*Nv))*exp((zeta+Ei_s_Ev_s)/(Kb*Tempr));
% tau_n1=(1/(sigma*Vth*Nc))*exp((Ec_s_Ei_s-zeta)/(Kb*Tempr));
% =========================================================================

Dit1=Dita+Ditd;
frequencies=logspace(3,6,100);
t1=5e-7;
t2=0.9/(Kb*Tempr);
zeta2=zeta;
tau_p1=t1*exp(t2*(zeta2));
tau_n1=t1*exp(-t2*(zeta2));
% tau_p1=1e3*ones(1,length(zeta));
% tau_n1=1e3*ones(1,length(zeta));

tau_p=interp1(zeta,tau_p1,Phi_s);
tau_n=interp1(zeta,tau_n1,Phi_s);
Dit=interp1(zeta,Dit1,Phi_s);

for ii = 1:length(frequencies)
    
    omega = 2*3.142*frequencies(ii);
    
    imw = 1i*omega;
    
    ap = 1./(imw*tau_p);
    an = 1./(imw*tau_n);
    
    A = sqrt( 1 + 2*(ap + an) + (ap - an).^2);
    B = -0.25*( 1 + 2*(ap+an) );
    
    Hg  = - (1./(imw*A)) .* log((1-2*A-4*B)./(1+2*A-4*B));
    
    Hcp = - (1/(2*imw)).*log(tau_n./tau_p) + 0.5*(1+ap-an).*Hg ;
    Hcn = - (1/(2*imw)).*log(tau_p./tau_n) + 0.5*(1+an-ap).*Hg ;
    
    Ggr = Q_SI*Dit.*(tau_n.*tau_p).^(-1).*Hg;
    
    Ctp = Q_SI*Dit.*tau_p.^(-1).*Hcp;
    Ctn = Q_SI*Dit.*tau_n.^(-1).*Hcn;
%     
    Ys(ii,:) = imw.*( HFCV.Cs' + Ctn + Ctp .* (Ggr ./ (imw.*Ctp+Ggr) ) );
    Gw(ii,:) = real(Ys(ii,:))./omega;
    
    Cap_Cs(ii,:) = imag(Ys(ii,:))./omega;
end

Cap=1./(1/Cox+1./Cap_Cs);


figure(1)
subplot(2,2,1)
semilogx(2*3.142*frequencies,Gw);
xlabel('\omega');
ylabel('G/\omega');
subplot(2,2,2)
semilogy(zeta*13.6,tau_p1,zeta*13.6,tau_n1)
legend('\tau_p','\tau_n');
xlabel('E-Ei');
ylabel('time constant (sec)');
subplot(2,2,3)
plot(Vg_new,Cap)
subplot(2,2,4)
semilogy(zeta*13.6,Ditd+Dita)
ylabel('Dit');
xlabel('E-Ei')


% =========================================================================

% After matching also plot Dit profile
% figure(2)
% semilogy(psi,Dita+Ditd,'sb');
% set(gca,'fontsize',15)
% xlabel('E-Ev (eV)');
% ylabel('D_{it} (cm^{-2} eV^{-1})');
% grid on;
% set(gca,'gridlinestyle','--');
% 
% % reply = input('please provide file name to write Dit profile \n reply "n" if you dont want exit\n','s');
% [fsave psave]=uiputfile('*.txt',Ideal_Dir);
% if ~(fsave==0)
%     fid=fopen([psave '/' fsave],'w');
%     fprintf(fid,'Vg \t C \n');
%     for i=1:length(Vg_new)
%         fprintf(fid,'%6.2f \t %e \n',Vg_new(i),C(i)*1e6);
%     end
%     fclose(fid);
% end



fprintf('See you\n');