function [rho, n_charge, p_charge, rho1] = charge_density(V_out,OutGenerate,Ec_new, Ev_new,OutSchr,mode)
% Function that calculates the chage densities based on either
% Schrodinger-Poisson or Poisson configurations 

global Kb hbar Q

    Ef          = OutGenerate.Ef;
    eta_C       = OutGenerate.eta_C;
    eta_V       = OutGenerate.eta_V;
    mn_eff      = OutGenerate.mn_eff; %For Spin Degeneracy
    mp_eff      = OutGenerate.mp_eff;
    Nd          = OutGenerate.Nd;
    Na          = OutGenerate.Na;
    Nc          = OutGenerate.Nc;
    Nv          = OutGenerate.Nv;
    Tempr       = OutGenerate.Tempr;
    No_of_nodes = length(OutGenerate.x);
    
    if mode == 1
        OutGenerate.control.Solvers.SP = 0;
    end


% if mode==0
%     if OutSchr.direction == 'l'
%     eta_C  = 2*ones(1,No_of_nodes);
%     mn_eff = OutGenerate.layer.layer_material.mn_t*ones(1,No_of_nodes);
%     elseif OutSchr.direction == 't'
%     eta_C = 4*ones(1,No_of_nodes);
%     mn_eff = OutGenerate.layer.layer_material.mn_l*ones(1,No_of_nodes);   
%     end
% end
    
    %persistent SR CLS IsFirst
    if(OutGenerate.control.Solvers.SP)
    Psi2_C      = OutSchr.Psi2_C;
    Psi2_V      = OutSchr.Psi2_V;
    Eigen_val_C = OutSchr.Eigen_val_C;
    Eigen_val_V = OutSchr.Eigen_val_V;

    SchrStart   = OutGenerate.SchrStart;
    SchrStop    = OutGenerate.SchrStop;
    UseSchr.C = 1; UseSchr.V = 1;

    else
       UseSchr.C = 0; UseSchr.V = 0; 
    end   
    

    n_charge = zeros(1,No_of_nodes);
    p_charge = zeros(1,No_of_nodes);
    
    IsFirst = 0;

    if(OutGenerate.control.Solvers.SP)
%     if isempty(IsFirst)
%         IsFirst = 0;
%     end

    if IsFirst == 0
        SR  = ~logical(1:No_of_nodes);
        SR(SchrStart:SchrStop) = true;
        CLS = ~SR;
    end

    IsFirst = IsFirst + 1;
    end

    %% n_charge

    if UseSchr.C
        
        temp = (max(Ef)-Eigen_val_C)/(Kb*Tempr);
        spin = 2*ones(1,No_of_nodes);
        % We are considering that many states, so that (Ef-E)/Kb*Tempr
        % of the highest state considered is less than ~ -8
        
        temp1 = find(temp<-8,1,'first');
        
        No_of_states_to_consider = 15;%temp1;
        
        for i=1:No_of_states_to_consider
            n_charge(SR) = n_charge(SR) + ( eta_C(SR).*mn_eff (SR).*Kb*Tempr/(pi*hbar^2) )...
                .* Psi2_C(:,i)' .* log( 1+exp( (Ef(SR)-Eigen_val_C(i) )/(Kb*Tempr) ) );
        end
        
        n_charge(CLS) = Nc(CLS).*exp((Ef(CLS)-Ec_new(CLS))/(Kb*Tempr) );
%         eta_FD=(Ef(CLS)-Ec_new(CLS))/(Kb*Tempr);
%         Nc_tmp=Nc(CLS);
%         zeta=linspace(0,10,1000);
%         n_charge_tmp=zeros(1,length(Nc_tmp));
%         for i=1:length(Ec_new(CLS))
%             n_charge_tmp(i)=Nc_tmp(i)*trapz(zeta,(zeta.^(0.5)./(1+exp(zeta-eta_FD(i)))));
%         end
%         n_charge(CLS)=n_charge_tmp;
    %         n_charge(CLS) = Nc(CLS).*(1./(1+exp(-(Ef(CLS)-Ec_new(CLS))/(Kb*Tempr))));
    end


    %% p_charge

    if  UseSchr.V
        
        temp = (Eigen_val_V-min(Ef))/(Kb*Tempr);
        
        % We are considering that many states, so that (E-Ef)/Kb*Tempr
        % of the highest state considered is less than ~ -8
        
        temp1 = find(temp<-8,1,'first');
        
        No_of_states_to_consider = 15;%temp1;
        
        for i=1:No_of_states_to_consider
            p_charge(SR) = p_charge(SR) + ( eta_V(SR).*mp_eff(SR).*Kb*Tempr/(pi*hbar^2) )...
                .* Psi2_V(:,i)' .* log( 1+exp( (Eigen_val_V(i)-Ef(SR) )/(Kb*Tempr) ) );
        end
        
        p_charge(CLS) = Nv(CLS).*exp((Ev_new(CLS)-Ef(CLS))/(Kb*Tempr) );
%         eta_FD=(Ev_new(CLS)-Ef(CLS))/(Kb*Tempr);
%         zeta=linspace(0,10,1000);
%         Nv_tmp=Nv(CLS);
%         p_charge_tmp=zeros(1,length(Nv_tmp));
%         for i=1:length(Ec_new(CLS))
%             p_charge_tmp(i)=Nv_tmp(i)*trapz(zeta,(zeta.^(0.5)./(1+exp(zeta-eta_FD(i)))));
%         end
%         p_charge(CLS)=p_charge_tmp;
    %      p_charge(CLS) = Nv(CLS).*(1./(1+exp(-(Ev_new(CLS)-Ef(CLS))/(Kb*Tempr))));
        
    end


% if mode==1
%     n_charge = Nc.*exp((Ef-Ec_new)/(Kb*Tempr)) ;
%     p_charge = Nv.*exp((Ev_new-Ef)/(Kb*Tempr) );
%     rho = Q*(Nd-Na+p_charge-n_charge);  % Q = |electron charge|
%     rho1= 0;
%     return;
% end

  
%%  Classical

if ~UseSchr.C
    % Classical electron Charge
    %n_charge = Nc.*exp((Ef-Ec_new)/(Kb*Tempr) );
%     figure(2);plot(Ef*13.6); hold all;plot(Ec_new*13.6);
%     rply=input();
%     'here'
      n_charge = Nc.*(1./(1+exp(-(Ef-Ec_new)/(Kb*Tempr))) );
%    eta_FD=(Ef-Ec_new)/(Kb*Tempr);    
%    zeta=linspace(0,10,1000);
%    for i=1:length(Ec_new)
%        n_charge(i)=Nc(i)*trapz(zeta,(zeta.^(0.5)./(1+exp(zeta-eta_FD(i)))));
%    end
end
if ~UseSchr.V
    % Classical hole charge
    %p_charge = Nv.*exp((Ev_new-Ef)/(Kb*Tempr) );
      p_charge = Nv.*(1./(1+exp(-(Ev_new-Ef)/(Kb*Tempr) )));
%    eta_FD=((Ev_new-Ef))/(Kb*Tempr);
%    zeta=linspace(0,10,1000);
%    for i=1:length(Ev_new)
%        p_charge(i)=Nv(i)*trapz(zeta,(zeta.^(0.5)./(1+exp(zeta-eta_FD(i)))));
%    end
end

rho = Q*(Nd-Na+p_charge-n_charge);  % Q = |electron charge|
rho1 = 0;




end
