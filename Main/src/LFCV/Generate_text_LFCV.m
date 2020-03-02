function Generate_text_LFCV(output_dir,file_prefix)

if ~exist(fullfile(output_dir,'mat_files','LFCV'),'dir')
    mkdir(fullfile(output_dir,'mat_files','LFCV'))
end

if output_dir(end)=='/'
    output_dir_LF=[output_dir(1:end-1) '_LF'];
else
    output_dir_LF=[output_dir '_LF'];
end

files = dir(fullfile(output_dir, 'mat_files','*.mat'));
M_CV = [];
for i=1:length(files)
   input_file = fullfile(output_dir, 'mat_files',files(i).name);
   load(input_file); 
       
   x  = OutGenerate.x;        
%    trap_loc_index = 11;
%    h1 = abs(x(trap_loc_index)-x(trap_loc_index-1));
%    h2 = abs(x(trap_loc_index+1)-x(trap_loc_index));
%    rho_it = Qit*2/(h1+h2);
   
%    rho(trap_loc_index) = rho(trap_loc_index) - rho_it;
   total_charge = (1.60217e-19 / sqrt(2) )*trapz(x,rho) * (0.529e-8)^-2;
%    total_charge = total_charge + Qit * (1.60217e-19 / sqrt(2) )*(0.529e-8)^-2;
   M_CV  = [M_CV; V_gate total_charge];
end

files = dir(fullfile(output_dir_LF, 'mat_files','*.mat'));
M_CV_LF = [];
for i=1:length(files)
   input_file = fullfile(output_dir_LF, 'mat_files',files(i).name);
   load(input_file); 
       
   x  = OutGenerate.x;        
%    trap_loc_index = 11;
%    h1 = abs(x(trap_loc_index)-x(trap_loc_index-1));
%    h2 = abs(x(trap_loc_index+1)-x(trap_loc_index));
%    rho_it = Qit*2/(h1+h2);
   
%    rho(trap_loc_index) = rho(trap_loc_index) - rho_it;
   total_charge = (1.60217e-19 / sqrt(2) )*trapz(x,rho) * (0.529e-8)^-2;
%    total_charge = total_charge + Qit * (1.60217e-19 / sqrt(2) )*(0.529e-8)^-2;
   M_CV_LF  = [M_CV_LF; V_gate total_charge];
end

M_CV = sortrows(M_CV,1);
M_CV_LF = sortrows(M_CV_LF,1);

dQ=M_CV(:,2)-M_CV_LF(:,2);
dV=M_CV(:,1)-M_CV_LF(:,1);

V=M_CV(:,1);
C = abs(dQ./dV);

filename = fullfile(output_dir, 'mat_files', 'LFCV', [file_prefix 'LFCV.mat']);
save(filename,'V','C');


end