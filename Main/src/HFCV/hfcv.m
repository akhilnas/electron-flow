% Newton - Raphson implementation of CVsimulator
% Started on 21 May 2014 by dhirendra

% Maintained by : Dhirendra (dhirendra22121987@gmail.com)

function hfcv(input_file)

OutGenerate=generate(input_file);


H=getH(OutGenerate);
H=sparse(H);


if ( OutGenerate.control.Solvers.SP || OutGenerate.control.Solvers.P )
    for i = 1:size(OutGenerate.surface.bias,2)
        try
            
            SPSCS_hf(H,OutGenerate,OutGenerate.surface.bias{i});
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
Generate_text_HFCV(output_dir,file_prefix);
generate_text_output_HF(output_dir);


end


% if OutGenerate.control.Solvers.Trap
%     
%     output_dir = OutGenerate.control.Output_directory;
%     files  = dir(fullfile(output_dir, 'mat_files','*.mat'));
%     
%     for i=1:length(files)
%         input_file = fullfile(output_dir, 'mat_files',files(i).name);
%         load(input_file);        
%     end
% end

fprintf('\n\nSimulation completed...\n');
end