function [Vg,Ctot,Gtot,Gpwtot,Cp,Gp,Gpw]=Dbt_admittance_p(HFCV,Simulation_INP)
% Global Constants
Kb = 1; Q_SI = 1.6e-19; Ab_CGS = 0.529e-8; Q=sqrt(2); m0=9.1e-31; hbar=1.054e-34;

% Dit admittance for p type substrate

f=Simulation_INP.f;
E=Simulation_INP.E;
x=Simulation_INP.x;
Dbt=Simulation_INP.Dbt;
tau_p0=Simulation_INP.tau_p0;
tau_n0=Simulation_INP.tau_n0;
Cox=Simulation_INP.Cox;
Ev_offset=Simulation_INP.Ev_offset;
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
% Ditd=Dit/2;
% Dita=Dit/2;
% % fprintf('setting acceptor and donor interface trap densities to be Dit/2\n');
% 
% for i=1:length(Vg)
%     fermi=1./(1+exp((E-Efermi(i))/Vt));
%     funct=Ditd.*(1-fermi)-Dita.*(fermi);
%     Qit(i)=Q_SI*trapz(E,funct);
% end
% delVg=-Qit'/Cox;
% Vg=Vg+delVg;


%% Admittance Calculation

% tau=0.5e-9;
% wt=2*pi*f*tau;
% funct=(wt.^2)./(1+wt.^2);
% Gt=1e6*funct;

tau=0.5e-9;
wt=2*pi*f*tau;
funct=1i*2*pi*f./(1+1i*wt);
Yt=1e6*funct;
% Ct=1e-3./(1+wt.^2);
% Gt_tmp=logspace(-6,6,100);
% ff=logspace(-1,10,100);
% Gt=interp1(ff,Gt_tmp,f);

ps=HFCV.ps;
ns=HFCV.ns;
alpha_p=HFCV.alpha_p;
Ctp_tmp=zeros(length(x),1);
for jj=1:length(f)
    imw=1i*Omega(jj);
    for ii=1:length(Vg)
        f0=1./(1+exp((E-Efermi(ii))/Vt));
        tau_p=tau_p0/(ps(ii));
%         tau_n=tau_n0/(ns(ii));
        itau_p=1./tau_p;
%         itau_n=1./tau_n;
         
       
        Evbar=Ev_offset*13.6+Efermi(ii);
        kappa=1e-2*sqrt(2*0.5*m0*Q_SI*Evbar)/hbar;
        
        for xi=1:length(x)
            
            Integral3=Dbt(xi,:).*itau_p*exp(-2*kappa*x(xi)).*f0.*(1-f0);
            Integral4=f0./(imw*f0.*(1-f0)+f0.*itau_p*exp(-2*kappa*x(xi)));
            Ctp_tmp(xi)=(Q_SI/Vt)*trapz(E,Integral3.*Integral4)*alpha_p(ii);
        end
        Ctp(ii)=trapz(x,Ctp_tmp);

    end
    
    Ys=imw*HFCV.Cs'+imw*Ctp;
    
    Cp(jj,:)=imag(Ys)/Omega(jj);
    Gp(jj,:)=imag(Ys)/Omega(jj);
    Gpw(jj,:)=real(Ys)/Omega(jj);
    
    Ytot=(1./(1/(imw*Cox+Yt(jj))+1./(Ys)));
    Ytot=1./(1./(Ytot));
    %(1/(Gs+imw*Cser)
    
  
    Ctot(jj,:)=imag(Ytot)/Omega(jj);
    Gpwtot(jj,:)=real(Ytot)/Omega(jj);
    Gtot(jj,:)=real(Ytot);
end

end