function [Dit, sigma_p]=Dit_sigma_gen(Fitting_Details,Dit_param,E,Vt)
if strcmp(Fitting_Details.Dit_fit_eq,'Gaussian') && strcmp(Fitting_Details.sigma_fit_eq,'Gaussian')
    Dit_params_ar=Dit_param(1:end-3*Fitting_Details.sigma_N_eq-1);
    sigma_params_ar=Dit_param(end-3*Fitting_Details.sigma_N_eq:end);
    
    Dit_params.Amp=Dit_params_ar(1:Fitting_Details.Dit_N_eq);
    Dit_params.m=Dit_params_ar(Fitting_Details.Dit_N_eq+1:2*Fitting_Details.Dit_N_eq);
    Dit_params.C=Dit_params_ar(2*Fitting_Details.Dit_N_eq+1:3*Fitting_Details.Dit_N_eq);
    Dit_params.base=Dit_params_ar(3*Fitting_Details.Dit_N_eq+1);
    
    sigma_params.Amp=sigma_params_ar(1:Fitting_Details.sigma_N_eq);
    sigma_params.m=sigma_params_ar(Fitting_Details.sigma_N_eq+1:2*Fitting_Details.sigma_N_eq);
    sigma_params.C=sigma_params_ar(2*Fitting_Details.sigma_N_eq+1:3*Fitting_Details.sigma_N_eq);
    sigma_params.base=sigma_params_ar(3*Fitting_Details.sigma_N_eq+1);
    
    Dit=0;
    sigma_p=0;
    
    for i=1:Fitting_Details.Dit_N_eq
        Dit=Dit+10^Dit_params.Amp(i)*exp(-(E-(Dit_params.m(i))).^2/(Dit_params.C(i)));
    end
    Dit=Dit+10^Dit_params.base;
    
    for i=1:Fitting_Details.sigma_N_eq
        sigma_p=sigma_p+10^sigma_params.Amp(i)*exp(-(E-(sigma_params.m(i))).^2/(sigma_params.C(i)));
    end
    sigma_p=sigma_p+10^sigma_params.base;
end

if strcmp(Fitting_Details.Dit_fit_eq,'Gaussian') && strcmp(Fitting_Details.sigma_fit_eq,'Exponential')
    Dit_params_ar=Dit_param(1:end-2*Fitting_Details.sigma_N_eq-1);
    sigma_params_ar=Dit_param(end-2*Fitting_Details.sigma_N_eq:end);
    
    Dit_params.Amp=Dit_params_ar(1:Fitting_Details.Dit_N_eq);
    Dit_params.m=Dit_params_ar(Fitting_Details.Dit_N_eq+1:2*Fitting_Details.Dit_N_eq);
    Dit_params.C=Dit_params_ar(2*Fitting_Details.Dit_N_eq+1:3*Fitting_Details.Dit_N_eq);
    Dit_params.base=Dit_params_ar(3*Fitting_Details.Dit_N_eq+1);
    
    sigma_params.Amp=sigma_params_ar(1:Fitting_Details.sigma_N_eq);
    sigma_params.alpha=sigma_params_ar(Fitting_Details.sigma_N_eq+1:2*Fitting_Details.sigma_N_eq);
    sigma_params.base=sigma_params_ar(2*Fitting_Details.sigma_N_eq+1);
    
    Dit=0;
    sigma_p=0;
    
    for i=1:Fitting_Details.Dit_N_eq
        Dit=Dit+10^Dit_params.Amp(i)*exp(-(E-(Dit_params.m(i))).^2/(Dit_params.C(i)));
    end
    Dit=Dit+10^Dit_params.base;
    
    for i=1:Fitting_Details.sigma_N_eq
        sigma_p=sigma_p+10^sigma_params.Amp(i)*exp(-sigma_params.alpha(i)*E/Vt);
    end
    sigma_p=sigma_p+10^sigma_params.base;
end

end