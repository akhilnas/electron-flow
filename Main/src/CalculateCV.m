%% CalculateCV
% Evaluates C-V of ideal C-V simulation and border trap simualtion

%%% function: CalculateCV
%  Input: OutGenerate, parsed input
% Output: none, fills output_dir
% Called by: main.m

function CalculateCV(OutGenerate)
% OG=OutGenerate;
% output_dir1='/pg/rs/vdhirendra/simulator_v2/device_lab/dit_scs/no_dit1';
%         files  = dir(fullfile(output_dir1, 'mat_files','*.mat'));
%         M_CV=[];
%         for i=1:length(files)
%             input_file = fullfile(output_dir1, 'mat_files',files(i).name);
%             load(input_file);
%             global Kb Q;
%             Tempr=OutGenerate.Tempr;
%             
%             Eg = 1.12/13.6;
%             Ei = Eg/2; % approximate intrinsic level
%             trap_loc_index = 101;
%             
%             zeta = linspace(0,Eg,50) - Ei; % Energy measured wrt to intrinsic level
%             
%             %Dita = zeros(1,length(zeta));
%             %Ditd = zeros(1,length(zeta));
%             
%             
%             DD = 5e12;
%             %Dita(round(0.5*length(zeta)):end) = DD;
%             %Ditd(1:round(0.5*length(zeta))-1) = DD;
%             
%             Dita = DD*ones(1,length(zeta));
%             Ditd = DD*ones(1,length(zeta));
%             
%             Ei_s = OutGenerate.Ei(trap_loc_index) - V_out(trap_loc_index);
%             Ef_s = OutGenerate.Ef(trap_loc_index);
%             
%             temp = 1 + exp( (zeta+Ei_s-Ef_s)/(Kb*Tempr) );
%             
%             Qit = Q*(0.529e-8)^2 * trapz(zeta*13.6,Ditd)  ...
%                 - Q*(0.529e-8)^2 * trapz(zeta*13.6,(Dita+Ditd)./temp);
%             
%             x = OutGenerate.x;
%             h1 = abs(x(trap_loc_index)-x(trap_loc_index-1));
%             h2 = abs(x(trap_loc_index+1)-x(trap_loc_index));
%             rho_it = Qit*2/(h1+h2);
%             %  rho
%             
%             
%             
%             M_CV=[M_CV; V_gate rho_it];
%         end
%         M_CV=sortrows(M_CV,1);
%         
%         
%         
%         plot(M_CV(:,2)); figure(2);
%         Qit_rho=M_CV(:,2);
% % Qit_rho=0;        
% 
% OutGenerate=OG;
if ( OutGenerate.control.Solvers.SP || OutGenerate.control.Solvers.P )
    for i = 1:size(OutGenerate.surface.bias,2)
        try
            
            SP_SCS(OutGenerate,OutGenerate.surface.bias{i});
        catch err
            idSegLast = regexp(err.identifier, '(?<=:)\w+$', 'match');
            if  strcmp(idSegLast,'ConvergenceFailed')
                fprintf('\n'); fprintf(err.message); fprintf('\n\n');
                return
            else
                rethrow(err);
            end
        end
    end

output_dir  = OutGenerate.control.Output_directory; 
file_prefix = OutGenerate.control.file_prefix;

fprintf('\nSaving the results in:\n %s \n',output_dir);
Generate_text_CV(output_dir,file_prefix);
generate_text_output(output_dir);

end


if OutGenerate.control.Solvers.Trap
    
    output_dir = OutGenerate.control.Output_directory;
    files  = dir(fullfile(output_dir, 'mat_files','*.mat'));
    
    for i=1:length(files)
        input_file = fullfile(output_dir, 'mat_files',files(i).name);
        load(input_file);        
    end
end

fprintf('\n\nSimulation completed...\n');
end