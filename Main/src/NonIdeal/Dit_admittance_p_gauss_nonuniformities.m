function [Vg,Ctot,Gtot,Gpwtot,Cp,Gp,Gpw]=Dit_admittance_p_gauss_nonuniformities(HFCV,Simulation_INP)
% Global Constants
Kb = 1; Q_SI = 1.6e-19; Ab_CGS = 0.529e-8; Q=sqrt(2);

% Dit admittance for p type substrate

f=Simulation_INP.f;
E=Simulation_INP.E;
Dit=Simulation_INP.Dit;
tau_p0=Simulation_INP.tau_p0;
tau_n0=Simulation_INP.tau_n0;
Cox=Simulation_INP.Cox;
Sigma_G=Simulation_INP.sigma_G;
sigma_G=Sigma_G*HFCV.Vt*13.6; % sigma for gaussian distribution

Efermi=HFCV.Efermi;
Vt=HFCV.Vt;


Omega     = 2*3.142*f;

%% Variable initialization
Vg=HFCV.Vg;
Ggr=zeros(1,length(Vg));
Ctn=zeros(1,length(Vg));
Ctp=zeros(1,length(Vg));
% Ys=zeros(length(f),length(Vg));
Ytot=zeros(length(f),length(Vg));
Ctot=zeros(length(f),length(Vg));
Gtot=zeros(length(f),length(Vg));
Gp=zeros(length(f),length(Vg));
Gpw=zeros(length(f),length(Vg));
Cp=zeros(length(f),length(Vg));
Gpwtot=zeros(length(f),length(Vg));
Qit=zeros(1,length(Vg));

%% =========== Qit addition ============
Ditd=Dit/2;
Dita=Dit/2;
% fprintf('setting acceptor and donor interface trap densities to be Dit/2\n');

for i=1:length(Vg)
    fermi=1./(1+exp((E-Efermi(i))/Vt));
    funct=Ditd.*(1-fermi)-Dita.*(fermi);
    Qit(i)=Q_SI*trapz(E,funct);
end
delVg=-Qit'/Cox;
Vg=Vg+delVg;

%% Admittance Calculation


Yit=zeros(1,10);




for jj=1:length(f)
    imw=1i*Omega(jj);
%     fprintf('jj = %d\n',jj);
    for ii=1:length(Vg)
        % _G stands for gaussian potential fluctuations 
        psi=linspace(HFCV.psi_s(ii)-3*sigma_G(ii),HFCV.psi_s(ii)+3*sigma_G(ii),10);
%         Cs_G=interp1(HFCV.psi_s,HFCV.Cs,psi,'linear','extrap');
        Efermi_G=interp1(HFCV.psi_s,HFCV.Efermi,psi,'linear','extrap');
%         ps_G=interp1(HFCV.psi_s,HFCV.ps,psi,'linear','extrap');
%         ns_G=interp1(HFCV.psi_s,HFCV.ns,psi,'linear','extrap');
%         alpha_p_G=interp1(HFCV.psi_s,HFCV.alpha_p,psi,'linear','extrap');
%         fprintf('ii = %d\n',ii);
        for kk=1:length(psi)
            
%             f0=1./(1+exp((E-Efermi_G(kk))/Vt));
%             tau_p=tau_p0/(ps_G(kk));
%             tau_n=tau_n0/(ns_G(kk));
%             itau_p=1./tau_p;
%             itau_n=1./tau_n;
%             
%             Integral1=Dit.*itau_p.*itau_n.*f0.*(1-f0);
%             Integral2=1./(imw*f0.*(1-f0)+f0.*itau_p+(1-f0).*itau_n);
%             
%             Integral3=Dit.*itau_p.*f0.*(1-f0);
%             Integral4=f0./(imw*f0.*(1-f0)+f0.*itau_p+(1-f0).*itau_n);
%             
%             Integral5=Dit.*itau_n.*f0.*(1-f0);
%             Integral6=(1-f0)./(imw*f0.*(1-f0)+f0.*itau_p+(1-f0).*itau_n);
%             
%             Ggr=(Q_SI/Vt)*trapz(E,Integral1.*Integral2)*alpha_p_G(kk);
%             Ctp=(Q_SI/Vt)*trapz(E,Integral3.*Integral4)*alpha_p_G(kk);
%             Ctn=(Q_SI/Vt)*trapz(E,Integral5.*Integral6)*alpha_p_G(kk);
%             
%             Y1=1/(1/(imw*Ctn)+1/(Ggr));
%             Yit(kk)=imw*Ctp+Y1;
        tau_p=interp1(E,tau_p0,Efermi_G(kk))./(ps);
%         tau_n=interp1(E,tau_n0,Efermi)./(ns);
%         itau_p=1./tau_p;
%         itau_n=1./tau_n;
        
        Integral5=(Vt)*interp1(E,Dit,Efermi_G(kk));
        Integral6=1./(imw*0.5*tau_p+1);

        Ctp=(Q_SI/Vt)*Integral5.*Integral6.*alpha_p;
        Yit(kk)=imw*Ctp;
        end
        gauss=normpdf(psi,HFCV.psi_s(ii),2*HFCV.Vt*13.6);
        Yit_avg=trapz(psi,gauss.*Yit);
        Ys=imw*HFCV.Cs(ii)+Yit_avg;
        Ytot=1./(1./(imw*Cox)+1./(Ys));
        Ctot(jj,ii)=imag(Ytot)/Omega(jj);
        Gtot(jj,ii)=real(Ytot);
        Gpwtot(jj,ii)=real(Ytot)/Omega(jj);
        
        Cp(jj,ii)=imag(Ys)/Omega(jj);
        Gp(jj,ii)=real(Ys);
        Gpw(jj,ii)=real(Ys)/Omega(jj);
    end
    
end

end