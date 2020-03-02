function [OutGenerate HFCV]=Ideal_HFCV(IdealFileDir, Tox)
Kb=1;

output_dir=IdealFileDir;
files  = dir(fullfile(output_dir, 'mat_files','*.mat'));

M_CV1=[];
for i=1:length(files)
    input_file = fullfile(output_dir, 'mat_files',files(i).name);
    load(input_file);
    
    x  = OutGenerate.x;
    [forget surface_loc_ind]=min(abs(OutGenerate.x-Tox));

    total_n_charge = -(1.60217e-19 )*trapz(x,n_charge) * (0.529e-8)^-2;
    total_p_charge =  (1.60217e-19 )*trapz(x,p_charge) * (0.529e-8)^-2;
    psi=(V0(surface_loc_ind))*19.2;
    PSI_S=(Ec_new(surface_loc_ind)-OutGenerate.Ec(surface_loc_ind))*13.6;
    ns=(n_charge(surface_loc_ind+0))* (0.529e-8)^-3; % surface electron density
    ps=(p_charge(surface_loc_ind+0))* (0.529e-8)^-3; % surface hole density
    
    
    M_CV1  = [M_CV1; V_gate psi total_p_charge total_n_charge ps ns];
end
    fprintf('surface location index is %d \n',surface_loc_ind);

M_CV1=sortrows(M_CV1,1);

output_dir=[output_dir '_HF'];
files  = dir(fullfile(output_dir, 'mat_files','*.mat'));
M_CV2=[];
for i=1:length(files)
    input_file = fullfile(output_dir, 'mat_files',files(i).name);
    load(input_file);
    
    x  = OutGenerate.x;

    total_n_charge = -(1.60217e-19 )*trapz(x,n_charge) * (0.529e-8)^-2;
    total_p_charge =  (1.60217e-19 )*trapz(x,p_charge) * (0.529e-8)^-2;
    psi=(V0(surface_loc_ind))*19.2;
    PSI_S=(Ec_new(surface_loc_ind)-OutGenerate.Ec(surface_loc_ind))*13.6;
    
    ns=(n_charge(surface_loc_ind+0))* (0.529e-8)^-3; % surface electron density
    ps=(p_charge(surface_loc_ind+0))* (0.529e-8)^-3; % surface hole density
    
    M_CV2  = [M_CV2; V_gate psi total_p_charge total_n_charge ps ns];
end

M_CV2=sortrows(M_CV2,1);

dVHF=M_CV2(:,1)-M_CV1(:,1);
dpsiHF=M_CV2(:,2)-M_CV1(:,2);
dpHF=M_CV2(:,3)-M_CV1(:,3);
dnHF=M_CV2(:,4)-M_CV1(:,4);

dQHF=-(dpHF+dnHF);

C_HF=(dQHF./dVHF);
Cs=(dQHF./dpsiHF);

Vg=M_CV1(:,1);
psi_s=M_CV1(:,2);

Qsemi=M_CV1(:,3);
Cox=-Qsemi(1)/(Vg(1)-psi_s(1));

ps=M_CV1(:,5);
ns=M_CV1(:,6);

% Calculation of alpha, the degeneracy correction factor
Vt=OutGenerate.Tempr*13.6;
dps=M_CV1(:,5)-M_CV2(:,5);
dns=M_CV1(:,6)-M_CV2(:,6);
alpha_p=(dps*Vt./ps)./dpsiHF;
alpha_n=-(dns*Vt./ns)./dpsiHF;


% Making Ev=0, all energies in eV
[forget,trap_loc_index]=min(abs(OutGenerate.x-Tox));
Ef_Ev = (OutGenerate.Ef(trap_loc_index)-OutGenerate.Ev(trap_loc_index))*13.6;  
% Ec = OutGenerate.Ec(trap_loc_index)*13.6+Ef_Ev;
% Ev = zeros(1,length(Ec));
% Eg = Ec-Ev;


Vt        = Kb*OutGenerate.Tempr*13.6; % thermal energy
Efermi    = psi_s+Ef_Ev;


HFCV.Vg=Vg;
HFCV.psi_s=psi_s;
HFCV.C=C_HF;
HFCV.Cs=Cs;
HFCV.Cox=Cox;
HFCV.ps=ps;
HFCV.ns=ns;
HFCV.alpha_p=alpha_p;
HFCV.alpha_n=alpha_n;
HFCV.Efermi=Efermi;
HFCV.Vt=Vt;

end