function Cox=getCox(H,OutGenerate)


% OutGenerate.surface.bias{1}(end+1)=OutGenerate.surface.bias{1}(end)+0.0001;
if exist('./tmp','dir');
    rmdir('./tmp/','s');
end
OutGenerate.control.Output_directory='./tmp/';
SPSCS_hf(H,OutGenerate,OutGenerate.surface.bias{1});


output_dir=OutGenerate.control.Output_directory;
if ~exist(fullfile(output_dir,'mat_files','CV'),'dir')
    mkdir(fullfile(output_dir,'mat_files','CV'))
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

M_CV = sortrows(M_CV,1);

dQ = diff(M_CV(:,2));
dV = diff(M_CV(:,1));

V = M_CV(2:end,1);
C = abs(dQ./dV);
Cox=max(C);
% rmdir('./tmp/','s');
end