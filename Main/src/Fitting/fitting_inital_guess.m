function [Dit_param0 ub lb]=fitting_inital_guess(Fitting_Details,E,Eg)
if strcmp(Fitting_Details.Dit_fit_eq,'Gaussian') && strcmp(Fitting_Details.sigma_fit_eq,'Gaussian')
    Amp=[];
    m=[];
    C=[];
    
    Dit_ini.Amp=log10(5e12)*ones(1,Fitting_Details.Dit_N_eq);
    M_tmp=linspace(min(E),max(E),Fitting_Details.Dit_N_eq+2);
    Dit_ini.m=M_tmp(2:end-1);
    Dit_ini.C=0.01*ones(1,Fitting_Details.Dit_N_eq);
    Dit_ini.base=log10(1e11);
    
    sigma_ini.Amp=log10(1e-16)*ones(1,Fitting_Details.sigma_N_eq);
    M_tmp=linspace(min(E),max(E),Fitting_Details.sigma_N_eq+2);
    sigma_ini.m=M_tmp(2:end-1);
    sigma_ini.C=0.01*ones(1,Fitting_Details.sigma_N_eq);
    sigma_ini.base=log10(1e-21);
    
    % ub
    ubd.Amp=log10(5e14)*ones(1,Fitting_Details.Dit_N_eq);
    ubd.m=(Eg+0.2)*ones(1,Fitting_Details.Dit_N_eq);
    ubd.C=100*ones(1,Fitting_Details.Dit_N_eq);
    ubd.base=log10(5e14);
    
    ubs.Amp=log10(1e-11)*ones(1,Fitting_Details.sigma_N_eq);
    ubs.m=(Eg+0.2)*ones(1,Fitting_Details.sigma_N_eq);
    ubs.C=100*ones(1,Fitting_Details.sigma_N_eq);
    ubs.base=log10(1e-11);
    
    %lb
    lbd.Amp=log10(1e10)*ones(1,Fitting_Details.Dit_N_eq);
    lbd.m=(-1)*ones(1,Fitting_Details.Dit_N_eq);
    lbd.C=0.00001*ones(1,Fitting_Details.Dit_N_eq);
    lbd.base=log10(1e11);
    
    lbs.Amp=log10(1e-21)*ones(1,Fitting_Details.sigma_N_eq);
    lbs.m=(-1)*ones(1,Fitting_Details.sigma_N_eq);
    lbs.C=0.00001*ones(1,Fitting_Details.sigma_N_eq);
    lbs.base=log10(1e-22);
    
    Dit_param0=[Dit_ini.Amp Dit_ini.m Dit_ini.C Dit_ini.base sigma_ini.Amp sigma_ini.m sigma_ini.C sigma_ini.base];
    ub=[ubd.Amp ubd.m ubd.C ubd.base ubs.Amp ubs.m ubs.C ubs.base];
    lb=[lbd.Amp lbd.m lbd.C lbd.base lbs.Amp lbs.m lbs.C lbs.base];
end

if strcmp(Fitting_Details.Dit_fit_eq,'Gaussian') && strcmp(Fitting_Details.sigma_fit_eq,'Exponential')
    Amp=[];
    m=[];
    C=[];
    
    Dit_ini.Amp=log10(5e12)*ones(1,Fitting_Details.Dit_N_eq);
    M_tmp=linspace(min(E),max(E),Fitting_Details.Dit_N_eq+2);
    Dit_ini.m=M_tmp(2:end-1);
    Dit_ini.C=0.01*ones(1,Fitting_Details.Dit_N_eq);
    Dit_ini.base=log10(1e11);
    
    sigma_ini.Amp=log10(1e-16)*ones(1,Fitting_Details.sigma_N_eq);
    alpha_tmp=linspace(min(E),max(E),Fitting_Details.sigma_N_eq+2);
    sigma_ini.alpha=alpha_tmp(2:end-1);
    sigma_ini.base=log10(1e-21);
    
    % ub
    ubd.Amp=log10(5e14)*ones(1,Fitting_Details.Dit_N_eq);
    ubd.m=(Eg+0.2)*ones(1,Fitting_Details.Dit_N_eq);
    ubd.C=100*ones(1,Fitting_Details.Dit_N_eq);
    ubd.base=log10(5e14);
    
    ubs.Amp=log10(1e-11)*ones(1,Fitting_Details.sigma_N_eq);
    ubs.alpha=(Eg+0.2)*ones(1,Fitting_Details.sigma_N_eq);
    ubs.base=log10(1e-11);
    
    %lb
    lbd.Amp=log10(1e10)*ones(1,Fitting_Details.Dit_N_eq);
    lbd.m=(-1)*ones(1,Fitting_Details.Dit_N_eq);
    lbd.C=0.00001*ones(1,Fitting_Details.Dit_N_eq);
    lbd.base=log10(5e11);
    
    lbs.Amp=log10(1e-21)*ones(1,Fitting_Details.sigma_N_eq);
    lbs.alpha=(-1)*ones(1,Fitting_Details.sigma_N_eq);
    lbs.base=log10(1e-21);
    
    Dit_param0=[Dit_ini.Amp Dit_ini.m Dit_ini.C Dit_ini.base sigma_ini.Amp sigma_ini.alpha sigma_ini.base];
    ub=[ubd.Amp ubd.m ubd.C ubd.base ubs.Amp ubs.alpha ubs.base];
    lb=[lbd.Amp lbd.m lbd.C lbd.base lbs.Amp lbs.alpha lbs.base];
end


end