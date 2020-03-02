function [OutGenerate LFCV]=Ideal_LFCV(IdealFileDir, Tox)
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
    M_CV1  = [M_CV1; V_gate psi total_p_charge total_n_charge];
end
    fprintf('surface location index is %d \n',surface_loc_ind);

M_CV1=sortrows(M_CV1,1);

output_dir=[output_dir '_LF'];
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
    M_CV2  = [M_CV2; V_gate psi total_p_charge total_n_charge];
end

M_CV2=sortrows(M_CV2,1);

dVLF=M_CV2(:,1)-M_CV1(:,1);
dpsiLF=M_CV2(:,2)-M_CV1(:,2);
dpLF=M_CV2(:,3)-M_CV1(:,3);
dnLF=M_CV2(:,4)-M_CV1(:,4);

dQLF=-(dpLF+dnLF);

C_LF=(dQLF./dVLF);
Cs=(dQLF./dpsiLF);

Vg=M_CV1(:,1);
psi_s=M_CV1(:,2);

Qsemi=M_CV1(:,3);

Cox=-Qsemi(1)/(Vg(1)-psi_s(1));

LFCV.Vg=Vg;
LFCV.psi_s=psi_s;
LFCV.C=C_LF;
LFCV.Cs=Cs;
LFCV.Cox=Cox;

end