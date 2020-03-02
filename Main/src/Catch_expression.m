function Catch_expression(err)
            idSegLast = regexp(err.identifier, '(?<=:)\w+$', 'match');
            
            switch idSegLast{:}
                case 'InvalidFileFid'
                    fprintf(err.message); fprintf('\n');
                case 'InvalidCommand'
                    fprintf(err.message); fprintf('\n');
                    %                     getReport(err,'extended');
                case 'FileNotFound'
                    fprintf(err.message); fprintf('\n');
                case 'InvalidValue'
                    fprintf(err.message); fprintf('\n');
                case 'OSError'
                    fprintf(err.message); fprintf('\n');
                case 'ConvergenceFailed'
                    fprintf(err.message); fprintf('\n');
                case 'ValueCheckError'
                    fprintf(err.message); fprintf('\n');
                case 'ValueNotSet'
                    fprintf(err.message); fprintf('\n');  
                case 'Unclassified'
                    fprintf(err.message); fprintf('\n');  
                otherwise
                    rethrow(err);
            end
end