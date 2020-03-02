function generate_text_output_LF(output_dir)


files = dir(fullfile(output_dir, 'mat_files','*.mat'));

% M_CV = zeros(length(files),2);

for i=1:length(files)
    input_file = fullfile(output_dir, 'mat_files',files(i).name);
    load(input_file);
    
    %-------------------------------------------------------------------------------------

    x  = OutGenerate.x;
    Nd = OutGenerate.Nd;
    Na = OutGenerate.Na;
    Ef = OutGenerate.Ef;
    if(OutGenerate.control.Solvers.SP)
        Psi2_C = OutSchr.Psi2_C;
    Psi2_V = OutSchr.Psi2_V;
    Eigen_val_V = OutSchr.Eigen_val_V;
    Eigen_val_C = OutSchr.Eigen_val_C;
    SchrStart   = OutGenerate.SchrStart;
    SchrStop    = OutGenerate.SchrStop;
    end
    
    %--------------------------------------------------------------------------------------
    [path, name, extension] = fileparts(input_file);
    
    out_dir = fullfile(output_dir, 'output_files');
    if ~exist(out_dir,'dir')
        mkdir(out_dir) ;
    end
    fid = fopen(fullfile(out_dir,[name '.txt']),'w');
    fprintf(fid,'Bias applied = %f Volts \n',V_gate);
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s \n','x(nm)','Ec(eV)','Ev(eV)','Ef(eV)','Nd(cm^-3)','Na(cm^-3)',...
        'n_charge(cm^-3)','p_charge(cm^-3)','rho(cm^-3) ' );
    fprintf(fid,'%.4f \t %.4f \t %.4f \t %.4f \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \n',...
        [0.0529*x; 13.6*Ec_new; 13.6*Ev_new; 13.6*Ef; (.529e-8)^-3 *Nd ; (.529e-8)^-3 *Na; ...
        (.529e-8)^-3 *n_charge;  (.529e-8)^-3 *p_charge; (.529e-8)^-3 *rho]);
    fclose(fid);
    
    wfn_dir = fullfile(output_dir, 'wfn_files');
    if ~exist(wfn_dir,'dir')
        mkdir(wfn_dir) ;
    end
   if(OutGenerate.control.Solvers.SP)
    
    fid = fopen(fullfile(wfn_dir, [name '_wfn' '.txt']),'w');
    fprintf(fid, '%s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \t %s \n ',...
        'x(nm)', 'Psi2_C_1','Psi2_C_2','Psi2_C_3','Psi2_C_4','Psi2_C_5','Psi2_C_6','Psi2_C_7','Psi2_C_8','Psi2_C_9','Psi2_C_10',...
        'Psi2_V_1','Psi2_V_2','Psi2_V_3','Psi2_V_4','Psi2_V_5','Psi2_V_6','Psi2_V_7','Psi2_V_8','Psi2_V_9','Psi2_V_10');
    
    temp = min(10,size(Psi2_C,1));
    fprintf(fid,['%.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e ' ...
        '\t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \t %.8e \n '],...
        [0.0529*x(SchrStart:SchrStop); Psi2_C(1:SchrStop-SchrStart+1,1:temp)';Psi2_V(1:SchrStop-SchrStart+1,1:temp)']);
    fprintf(fid,'\n\n\n\n');
    fprintf(fid,'%s \t %s \n','Cond Band Eigen Val(eV)','valence Band Eigen Val(eV)');
    
    temp = min(30,length(Eigen_val_C));
    fprintf(fid,'%.8f \t \t \t \t \t %.8f \n',[13.6*Eigen_val_C(1:temp)';13.6*Eigen_val_V(1:temp)']);
    
    fclose(fid);
end
    
%     % CV
%     total_charge = (1.60217e-19 / sqrt(2) )*trapz(x,rho) * (0.529e-8)^-2;
%     M_CV(i,:) = [V_gate total_charge];
end

% M_CV = sortrows(M_CV,1);
% 
% dQ = diff(M_CV(:,2));
% dV = diff(M_CV(:,1));
% C  = abs(dQ./dV);

files = dir(fullfile(output_dir, 'mat_files','LFCV','*.mat'));
input_file = fullfile(output_dir, 'mat_files','LFCV',files.name);
load(input_file);

CV_dir = fullfile(output_dir, 'CV_files');
if ~exist(CV_dir,'dir')
    mkdir(CV_dir) ;
end
fid = fopen(fullfile(CV_dir, 'LFCV.txt'),'w');
fprintf(fid,'%s \t %s \n','V (Volts)','C (mic F cm^-2)');
fprintf(fid,'%.4f \t %.8f \n',[V';1e6*C']);
fclose(fid);

end
