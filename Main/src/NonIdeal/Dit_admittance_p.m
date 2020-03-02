function [Vg,Ctot,Gtot,Gpwtot,Cp,Gp,Gpw]=Dit_admittance_p(HFCV,Simulation_INP)
% Global Constants
Kb = 1; Q_SI = 1.6e-19; Ab_CGS = 0.529e-8; Q=sqrt(2);

% Dit admittance for p type substrate

f=Simulation_INP.f;
E=Simulation_INP.E;
Dit=Simulation_INP.Dit;
tau_p0=Simulation_INP.tau_p0;
tau_n0=Simulation_INP.tau_n0;
Cox=Simulation_INP.Cox;
% Gs=Simulation_INP.Gs;
% Cser=Simulation_INP.Cser;

Efermi=HFCV.Efermi;
Vt=HFCV.Vt;


Omega     = 2*3.142*f;

%% Variable initialization
Vg=HFCV.Vg;
Ggr=zeros(1,length(Vg));
Ctn=zeros(1,length(Vg));
Ctp=zeros(1,length(Vg));
Ys=zeros(length(f),length(Vg));
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

ps=HFCV.ps;
ns=HFCV.ns;
alpha_p=ones(1,length(Vg));%

for jj=1:length(f)
    imw=1i*Omega(jj);
    for ii=1:length(Vg)
        f0=1./(1+exp((E-Efermi(ii))/Vt));
        tau_p=tau_p0/(ps(ii));
        tau_n=tau_n0/(ns(ii));
        itau_p=1./tau_p;
        itau_n=1./tau_n;
        
        Integral1=Dit.*itau_p.*itau_n.*f0.*(1-f0);
        Integral2=1./(imw*f0.*(1-f0)+f0.*itau_p+(1-f0).*itau_n);
        
        Integral3=Dit.*itau_p.*f0.*(1-f0);
        Integral4=f0./(imw*f0.*(1-f0)+f0.*itau_p+(1-f0).*itau_n);
        
        Integral5=Dit.*itau_n.*f0.*(1-f0);
        Integral6=(1-f0)./(imw*f0.*(1-f0)+f0.*itau_p+(1-f0).*itau_n);
        
        Ggr(ii)=(Q_SI/Vt)*trapz(E,Integral1.*Integral2)*alpha_p(ii);
        Ctp(ii)=(Q_SI/Vt)*trapz(E,Integral3.*Integral4)*alpha_p(ii);
        Ctn(ii)=(Q_SI/Vt)*trapz(E,Integral5.*Integral6)*alpha_p(ii);

        
    end
    
    
    Y1=1./(1./(imw*Ctn)+1./(Ggr));

    Ys=imw*HFCV.Cs'+imw*Ctp;%+Y1;
    
    Cp(jj,:)=imag(Ys)/Omega(jj);
    Gp(jj,:)=imag(Ys)/Omega(jj);
    Gpw(jj,:)=real(Ys)/Omega(jj);
    
        Cit(jj,:)=imag(Ys)/Omega(jj);
    Git(jj,:)=imag(Ys)/Omega(jj);
    Gpwit(jj,:)=real(Ys)/Omega(jj);
    
    Ytot=(1./(1/(imw*Cox)+1./(Ys)));
    Ytot=1./(1./(Ytot));
    %(1/(Gs+imw*Cser)
    
  
    Ctot(jj,:)=imag(Ytot)/Omega(jj);
    Gpwtot(jj,:)=real(Ytot)/Omega(jj);
    Gtot(jj,:)=real(Ytot);
end

end