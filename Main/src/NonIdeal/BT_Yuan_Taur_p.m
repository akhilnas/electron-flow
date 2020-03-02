function [Vg,Ctot,Gtot,Gpwtot,Cp,Gp,Gpw]=BT_Yuan_Taur_p(HFCV,Simulation_INP)
% Global Constants
Kb = 1; Q_SI = 1.6e-19; Ab_CGS = 0.529e-8; Q=sqrt(2); m0=9.1e-31; hbar=1.054e-34;

% Dit admittance for p type substrate

f=Simulation_INP.f;
E=Simulation_INP.E;
x=Simulation_INP.x;
eps_ox=Simulation_INP.eps_ox;
Dbt=Simulation_INP.Dbt;
tau_p0=Simulation_INP.tau_p0;
tau_n0=Simulation_INP.tau_n0;
Cox=Simulation_INP.Cox;
Ev_offset=Simulation_INP.Ev_offset;
Ec_offset=Simulation_INP.Ec_offset;
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

ps=HFCV.ps;
% ps=HFCV.ps;
% alpha_p=HFCV.alpha_p;


% Ctpbt=zeros(1,length(Vg));
% 
Dbt1=zeros(length(x),length(Vg));
for xi=1:length(x)
    Dbt1(xi,:)=interp1(E,Dbt(xi,:),Efermi);
end
tau_n1=interp1(E,tau_n0,Efermi)./ps;
% 
% 
% for jj=1:length(f)
%     imw=1i*Omega(jj);
%     
%     for ii=1:length(Vg)
%         Evbar=Ev_offset*13.6+Efermi(ii);
%         kappa=1e-2*sqrt(2*0.5*m0*Q_SI*Evbar)/hbar;
%         Ctpbt(ii)=alpha_p(ii)*Q_SI*trapz(x,Dbt1(:,ii)'./(1+0.5*imw*tau_p1(ii)*exp(2*kappa*x)));
%     end
%     
%     Ys=imw*HFCV.Cs'+imw*Ctpbt;
%     
%     Cp(jj,:)=imag(Ys)/Omega(jj);
%     Gp(jj,:)=imag(Ys)/Omega(jj);
%     Gpw(jj,:)=real(Ys)/Omega(jj);
%     
%     Ytot=(1./(1/(imw*Cox)+1./(Ys)));
%     Ytot=1./(1./(Ytot));
%     
%   
%     Ctot(jj,:)=imag(Ytot)/Omega(jj);
%     Gpwtot(jj,:)=real(Ytot)/Omega(jj);
%     Gtot(jj,:)=real(Ytot);
% end

for Vgind=1:length(Vg)

        Evbar=Ev_offset*13.6+Efermi(Vgind);
        kappa=1e-2*sqrt(2*0.5*m0*Q_SI*Evbar)/hbar;
    
    for jj=1:length(f)
        omega=2*pi*f(jj);
        imw=1i*omega;
        Y=zeros(1,length(x));
        Y(1)=imw*HFCV.Cs(Vgind);
        for ii=1:length(x)-1
            %      Y(ii+1)=(-Y(ii)*Y(ii)/(imw*eps_ox)+imw*Q_SI*Nbt(ii)/(1+imw*tau_0*exp(2*kappa*x(ii))))*(x(ii+1)-x(ii))+Y(ii);
            %  Y(ii+1)=(-Y(ii)*Y(ii)/(imw*eps_ox)+imw*Q_SI*Nbt(ii)/(1+imw*tau_0*exp(2*kappa*x(ii))))*(x(ii+1)-x(ii))+Y(ii);
%             Y(ii+1)=(-Y(ii)*Y(ii)/(imw*eps_ox(xi))+Q_SI*Dbt(ii,Vgind)'*log(1+imw*tau_n1(Vgind)*exp(2*kappa*x(ii)))/(tau_n1(Vgind)*exp(2*kappa*x(ii))))*(x(ii+1)-x(ii))+Y(ii); %2012 paper model
            Y(ii+1)=Y(ii)-Y(ii)*Y(ii)*(x(ii+1)-x(ii))/(imw*eps_ox(ii))+(Q_SI*Dbt(ii,Vgind)*log(1+imw*tau_n1(Vgind)*exp(2*kappa*x(ii)))/(tau_n1(Vgind)*exp(2*kappa*x(ii))))*(x(ii+1)-x(ii));
            % Y(ii+1)=(-Y(ii)*Y(ii)/(imw*eps_ox))*(x(ii+1)-x(ii))+Y(ii);
        end
        Ytot(jj,Vgind)=Y(end);
        Ctot(jj,Vgind)=imag(Y(end)/omega);
        Gpwtot(jj,Vgind)=real(Y(end)/omega);
        Gtot(jj,Vgind)=real(Y(end));
    end

end


end