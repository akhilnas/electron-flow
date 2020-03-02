%% MAIN FUNCTION
%  Displays initial text on terminal, sets default_path by accessing 1DPS_pref.txt
% Displays options for simulation, calls major files from this
% functions called

function  main()
addpath(genpath('/pg/rs/vdhirendra/simulator_v3/src/'))

fprintf('\n\n');
fprintf('       1D SCHRODINGER POISSON SOLVER           \n\n')
fprintf('    Department of Electrical Engineering         \n')
fprintf('                IIT Bombay                       \n')
fprintf('               Powai, Mumbai                     \n')
fprintf('                 400 076                         \n')
fprintf('---------------------------------------------\n\n\n')

if isdeployed()
    warning('off', 'MATLAB:nearlySingularMatrix');
end

%% Setting default_path



while(1)
    
    try
        WORK_DIR = parser_default_path()           
        break ;
    catch err
        idSegLast = regexp(err.identifier, '(?<=:)\w+$', 'match');
        
        if  (strcmp(idSegLast, 'InvalidFileFid') || strcmp(idSegLast, 'FolderNotExist') || strcmp(idSegLast, 'InvalidCommand'))
            fprintf(err.message); fprintf('\n');
            while(1)
                fprintf('(1) Please rectify the error and hit (c) and (Enter) to continue \n');
                fprintf('(2) Please hit (q) and (Enter) to exit \n\n');
                reply = input(' ','s');
                if strcmp(reply,'c') || strcmp(reply,'C')
                    break;
                elseif strcmp(reply,'q') || strcmp(reply,'Q')
                    fprintf('Program exiting.. \nSee you again.. Have a good day.. bye.. :) \n\n');
                    return;
                elseif strcmp(reply,'clc')
                    clc;
                end
            end
        end
    end
    
end

%% Displays Simulation options

while(1)
    
    fprintf('\n');
    fprintf('(1) Type (1) to run ideal CV simulation \n');
    fprintf('(2) Type (2) to run freq dependent inversion simulation \n');
    fprintf('(3) Type (3) to run Interface Traps CV simulation \n');
    fprintf('(4) Type (4) to run Border Traps CV simulation \n');
    fprintf('(5) Type (q) to terminate the program \n');
    fprintf('(6) Type (6) to plot ideal C-V \n');
    fprintf('Type valid function expression you want to evaluate in custom \n');
    
    reply = input(' ','s');
    if isempty(reply)OutGenerate=generate(WORK_DIR);
        continue;
    elseif strcmp(reply,'clc')
        clc;
    elseif strcmp(reply,'q') || strcmp(reply,'Q')
        fprintf('Program exiting.. \nSee you again.. Have a good day.. bye.. :) \n\n');
        break;
%%  If ideal C-V is selected, calles 'generate' to parse the input file		
    elseif strcmp(reply,'1') 
        
        try
            clear functions ;
            
            OutGenerate = generate(WORK_DIR);
        catch err
            Catch_expression(err);
            continue;
        end
        
        fprintf('Starting the computation...          \n\n\n\n');
 %% Uses 'Calculate CV' to generate ideal C-V       
        tic
        CalculateCV(OutGenerate);
        toc
        
        fprintf('\n\n\n');
        
    elseif strcmp(reply,'3')
        
         while(1)
             fprintf('(1) Type (1) to give Dit input file\n');
             fprintf('(2) Type (2) to input constant profile \n');
             fprintf('(3) Type (3) to  input Gaussian profile\n');  
             fprintf('(4) Type (q) to  go to main interface');  
            reply = input(' ','s');
            if isempty(reply)
              continue;
               elseif strcmp(reply,'clc')
               clc;
               elseif strcmp(reply,'q') || strcmp(reply,'Q')
                fprintf('out of this');
                break;
            elseif strcmp(reply,'1') 
                try
            clear functions;
  %% dit simulation: 'Calculate_CV_Dit' (add html )                    
            Calculate_CV_Dit(WORK_DIR);
            fprintf('Simulation completed successfully!');
            fprintf('\n\n\n');
            catch err
            Catch_expression(err);
            continue;
                end
            elseif strcmp(reply,'2') 
                try
            clear functions;
  %% dit simulation: 'Calculate_CV_Dit' (add html )                    
            Calculate_CV_Dit_constant(WORK_DIR);
            fprintf('Simulation completed successfully!');
            fprintf('\n\n\n');
            catch err
            Catch_expression(err);
            continue;
                end
            elseif strcmp(reply,'3') 
                try
            clear functions;
  %% dit simulation: 'Calculate_CV_Dit' (add html )                    
            Calculate_CV_Dit_gauss(WORK_DIR);
            fprintf('Simulation completed successfully!');
            fprintf('\n\n\n');
            catch err
            Catch_expression(err);
            continue;
                end    
         else
        fprintf('Please provide a valid input \n');
    
            end
            end
        
        
    elseif strcmp(reply,'4')
        
        try
            clear functions;
     %% BT simulation 'Calculate_CV_BT' (add html)                    
            Calculate_CV_BT(WORK_DIR);
            
            fprintf('Simulation completed successfully!');
            fprintf('\n\n\n');
        catch err
            Catch_expression(err);
            continue;
        end
        
     elseif strcmp(reply,'2')
        
        try
            clear functions;
                        
            Calculate_CV_Freq(WORK_DIR);
            
            fprintf('Simulation completed successfully!');
            fprintf('\n\n\n');
        catch err
            Catch_expression(err);
            continue;
        end
        
    elseif strcmp(reply,'6')
        
        try
            clear functions;
            if(exist('OutGenerate'))
                fprintf('Plotting ideal CV!');
                fprintf('\n\n\n');
                plotCV(OutGenerate.control.Output_directory)
            else
                fprintf('Please first hit enter and read controls from file'); 
                fprintf('\n\n\n');
            end
        catch err
            Catch_expression(err);
            continue;
        end
        
    else
        try
        eval(reply)
        catch err
            err
        fprintf('Please provide a valid input \n');
        end
    end
    
    
end 


end 

