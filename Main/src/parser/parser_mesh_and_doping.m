%% TODO: where and what
            j = 1;
            while (~strcmp(C{1}{i+j},'end'))
             switch C{1}{i+j}
                    case 'layer'
                        
                            for p=1:size(var_layer,2)
                                if (strcmp(eval([var_layer{p} '.layer_name']) , C{1}{i+j+1} ) )
                                    temp_layer =  var_layer{p};
                                    break;
                                end
                            end
                            
                        k = 2;
                        
                        sublayer_no = 0; 
                        var_sublayer = {};
                        while (~strcmp(C{1}{i+j+k},'end'))

                            switch C{1}{i+j+k}
                                case 'sublayer'
                                    var_sublayer = {var_sublayer{:} genvarname('sublayer',var_sublayer)};
                                    eval([var_sublayer{end} '= ClassSubLayer;']);             
                                    sublayer_no = sublayer_no +1;             
                                    eval([var_sublayer{sublayer_no} '.sublayer_name =' char(39) C{1}{i+j+k+1} char(39) ';']);
             
                                    l = 2;
                                    while (~strcmp(C{1}{i+j+k+l},'end'))
                                        switch C{1}{i+j+k+l}
                                            case 'sublayer_begin'
                                                eval([var_sublayer{sublayer_no} '.Pos_begin =' C{1}{i+j+k+l+1} ';']);
                                            case 'sublayer_end'
                                                eval([var_sublayer{sublayer_no} '.Pos_end ='   C{1}{i+j+k+l+1} ';']);
                                            case 'grid_spacing'
                                                eval([var_sublayer{sublayer_no} '.grid_spacing =' C{1}{i+j+k+l+1} ';']);
                                            case 'Nd'
                                                eval([var_sublayer{sublayer_no} '.Nd =' C{1}{i+j+k+l+1} ';']);
                                            case 'Na'
                                                eval([var_sublayer{sublayer_no} '.Na =' C{1}{i+j+k+l+1} ';']);
                                            case '<#'                                                
                                            % code for comment skipping
                                                m = 1;
                                                while ~strcmp(C{1}{i+j+k+l+m},'#>')
                                                    m = m+1;
                                                end
                                                l = l+m;
                                            otherwise
                                                %error
                                        end
                                        
                                        l = l+1;
                                    end
                                    k = k+l;
                                case 'grid_spacing'
                                    eval([temp_layer '.grid_spacing =' C{1}{i+j+k+1} ';']); 
                                case 'Na'
                                    eval([temp_layer '.Na =' C{1}{i+j+k+1} ';']); 
                                case 'Nd'
                                    eval([temp_layer '.Nd =' C{1}{i+j+k+1} ';']); 
                                case '<#'
                                    m = 1;
                                    while ~strcmp(C{1}{i+j+k+m},'#>')
                                        m = m+1;
                                    end
                                    k = k+m;                                    
                                otherwise
                                    %error
                            end                            
                            
                            k = k+1;
                        end
                    AA =[];
                    for pp=1:size(var_sublayer,2)
                        AA = [AA eval(var_sublayer{pp})];                         
                    end
                    eval([temp_layer '.sub_layers = AA;']);
                    
                    j = j+k;
                    %----------------------------------------------------------------
                case '<#'
                   %comment
                    m = 1;
                    while ~strcmp(C{1}{i+j+m},'#>')
                         m = m+1;
                    end
                    j = j+m;
                 otherwise 
                     %error
             end
             j = j+1;
            end
            i = i+j;
        %----------------------------------------------------------------------------------