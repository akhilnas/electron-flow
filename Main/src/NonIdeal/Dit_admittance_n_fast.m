function [Vg,Ctot,Gtot,Gpwtot,Cp,Gp,Gpw]=Dit_admittance_n_fast(HFCV,Simulation_INP)
% Global Constants
Kb = 1; Q_SI = 1.6e-19; Ab_CGS = 0.529e-8; Q=sqrt(2);

% Dit admittance for p type substrate

f=Simulation_INP.f;
E=Simulation_INP.E;
Dit=Simulation_INP.Dit;
tau_n0=Simulation_INP.tau_n0;
Cox=Simulation_INP.Cox;
% Gs=Simulation_INP.Gs;
% Cser=Simulation_INP.Cser;

Efermi=HFCV.Efermi;
Vt=HFCV.Vt;


Omega     = 2*3.142*f;

%% Variable initialization
Vg=HFCV.Vg;
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
alpha_n=HFCV.alpha_n;
ns=HFCV.ns;
for jj=1:length(f)
    imw=1i*Omega(jj);


%         tau_p=interp1(E,tau_p0,Efermi)./(ns);
        tau_n=interp1(E,tau_n0,Efermi)./(ns);
%         itau_p=1./tau_p;
%         itau_n=1./tau_n;
        
        Integral5=(Vt)*interp1(E,Dit,Efermi);
        Integral6=1./(imw*0.5*tau_n+1);

        Ctn=(Q_SI/Vt)*Integral5.*Integral6.*alpha_n;

        
    Ys=imw*HFCV.Cs+imw*Ctn;
    
    Cp(jj,:)=imag(Ys)/Omega(jj);
    Gp(jj,:)=imag(Ys)/Omega(jj);
    Gpw(jj,:)=real(Ys)/Omega(jj);
    
    Ytot=(1./(1/(imw*Cox)+1./(Ys)));
    Ytot=1./(1./(Ytot));
    %(1/(Gs+imw*Cser)
    
  
    Ctot(jj,:)=imag(Ytot)/Omega(jj);
    Gpwtot(jj,:)=real(Ytot)/Omega(jj);
    Gtot(jj,:)=real(Ytot);
end

end