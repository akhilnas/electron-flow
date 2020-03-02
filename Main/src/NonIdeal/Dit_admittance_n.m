function [Vg,Ctot,Gtot,Gpwtot,Cp,Gp,Gpw]=Dit_admittance_n(HFCV,Simulation_INP)
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
Ctn=zeros(1,length(Vg));
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

tau=1e-9;
wt=2*pi*f*tau;
funct=(wt.^2)./(1+wt.^2);
Gt=0.5e6*funct;
% Gt_tmp=logspace(-6,6,100);
% ff=logspace(-1,10,100);
% Gt=interp1(ff,Gt_tmp,f);

Gs=(1/0.01)/3600e-8;


ps=HFCV.ps;
ns=HFCV.ns;
alpha_n=HFCV.alpha_n;

for jj=1:length(f)
    imw=1i*Omega(jj);
    for ii=1:length(Vg)
        f0=1./(1+exp((E-Efermi(ii))/Vt));
        tau_p=tau_p0/(ps(ii));
        tau_n=tau_n0/(ns(ii));
        itau_p=1./tau_p;
        
        itau_n=1./tau_n;
        Integral5=Dit.*itau_n.*f0.*(1-f0);
        Integral6=(1-f0)./(imw*f0.*(1-f0)+f0.*itau_p+(1-f0).*itau_n);

        Ea=E(~isnan(Integral6));
        Integral5a=Integral5(~isnan(Integral6));
        Integral6a=Integral6(~isnan(Integral6));
        

        
        Ctn(ii)=(Q_SI/Vt)*trapz(Ea,Integral5a.*Integral6a)*alpha_n(ii);
                
    end

    Ys=imw*HFCV.Cs'+imw*Ctn;
    
    Cp(jj,:)=imag(Ys)/Omega(jj);
    Gp(jj,:)=imag(Ys)/Omega(jj);
    Gpw(jj,:)=real(Ys)/Omega(jj);
    
    Ytot=(1./(1/(imw*Cox+Gt(jj))+1./(Ys)));
    Ytot=1./(1./(Ytot));
  
    Ctot(jj,:)=imag(Ytot)/Omega(jj);
    Gpwtot(jj,:)=real(Ytot)/Omega(jj);
    Gtot(jj,:)=real(Ytot);
end

end