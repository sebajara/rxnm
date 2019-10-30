function [derscript,indS] = model2der(model,STR1,STR1val,STR2,STR2val)
% Takes a MODEL struct and creates a Matlab derivative function
% suitable for matlab ode solvers. It return the name of the created
% script.
%
% [DERSCRIPT] = model2der(MODEL)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
    OPTIONAL_STRS = {'Mode','Compile'};
    mode = 0; % 0 = matlab, 1 = ODEMEX
    compile = 0;
    if nargin < 1
        error('Model must be specified');
    elseif ((nargin == 2) | (nargin == 4) | (nargin == 6) | (nargin == 8))
        error('Optional inputs must be accompanied by values');
    elseif nargin == 3
        Input_Values = {STR1, STR1val};
    elseif nargin == 5
        Input_Values = {STR1, STR1val;STR2, STR2val};
    elseif nargin > 5
        error('At most 5 inputs are allowed');
    end
    for i=1:1:((nargin - 1)/2)
        STR = Input_Values{i,1};
        STRval = Input_Values{i,2};
        z = strmatch(STR,OPTIONAL_STRS,'exact');
        if z == 1
            mode = STRval;
        elseif z == 2
            compile = STRval;
        else
            error ('WTH!');
        end
    end
    
    rxns = fieldnames(model.rxns);
    species = fieldnames(model.species);
    nrxns = length(rxns);
    nspec = length(species);
    npara = length(fieldnames(model.params));
    nvars = length(fieldnames(model.vars));
    SM = stoichmatrix(model);
    
    indS = struct;
    indS.p = struct;
    indS.s = struct;
    indS.u = struct;
    indS.c = struct;
    
    if(mode==0) %% Normal MATLAB derivative file, suitable for ode15s
        derscript = [model.name,'_MATLABder'];
        
        stream = fopen([derscript,'.m'],'w');

        fprintf(stream,'function [dxdt] = derivative (t,x,p)\n');
        fprintf(stream,'%% Derivative for the reaction model %s\n',model.name);
        
        fprintf(stream,'\n');
        for i = 1:npara
            fprintf(stream,'   %s = p(%s);\n',model.ind2par{i},num2str(i));
        end
        fprintf(stream,'\n');
        for i = 1:nspec
            fprintf(stream,'   %s = x(%s);\n',model.ind2spc{i},num2str(i));
        end
        fprintf(stream,'\n');
        for i = 1:nvars
            formula = model.vars.(model.ind2var{i}).formula;
            fprintf(stream, '   %s = %s;\n',model.ind2var{i},formula);
        end
        fprintf(stream,'\n');
        
        fprintf(stream,'\n');
        fprintf(stream,'   dxdt = zeros(size(x));\n');
        
        for j = 1:nspec
            fprintf(stream, '   dxdt(%d) = 0',j);
            for i = 1:nrxns
                if (~(SM(i,j)==0))
                    rate = model.rxns.(rxns{i}).rate;
                    c = num2str(SM(i,j));
                    if (SM(i,j)>0)
                        c = ['+',c];
                    end
                    fprintf(stream, ' %s*(%s)',c,rate);
                end
            end
            fprintf(stream,[';%%   %s\n'],species{j}); 
        end
        fprintf(stream,'end\n');
        fprintf(stream,'\n');
    else %% Suitable for ODEMEX, and will be compiled later
        derscript = [model.name,'_ODEMEXder'];
        
        stream = fopen([derscript,'.m'],'w');
                
        fprintf(stream,'function [dxdt] = derivative (t,x,p,u,indS)\n');
        fprintf(stream,'%% Derivative for the reaction model ODEMEX%s\n',model.name);
        
        fprintf(stream,'\n');
        for i = 1:npara
            fprintf(stream,'   %s = p(indS.p.%s);\n',model.ind2par{i},model.ind2par{i});
            indS.p.(model.ind2par{i}) = i;
        end
        fprintf(stream,'\n');
        for i = 1:nspec
            fprintf(stream,'   %s = x(indS.s.%s);\n',model.ind2spc{i},model.ind2spc{i});
            indS.s.(model.ind2spc{i}) = i;
        end
        fprintf(stream,'\n');
        for i = 1:nvars
            formula = model.vars.(model.ind2var{i}).formula;
            fprintf(stream, '   %s = %s;\n',model.ind2var{i},formula);
        end
        fprintf(stream,'\n');
        
        fprintf(stream,'\n');
        %fprintf(stream,'   dxdt = zeros(%d,1);\n',nspec);
        
        for j = 1:nspec
            fprintf(stream, '   dxdt(indS.s.%s) = 0',model.ind2spc{j});
            for i = 1:nrxns
                if (~(SM(i,j)==0))
                    rate = model.rxns.(rxns{i}).rate;
                    if(regmbool(rate,'\^')) 
                        % replace ^ by pow
                        % otherwise won't compile 
                        rate = odemexsynrate(rate);
                    end
                    c = num2str(SM(i,j));
                    if (SM(i,j)>0)
                        c = ['+',c];
                    end
                    fprintf(stream, ' %s*(%s)',c,rate);
                end
            end
            fprintf(stream,[';%%   %s\n'],species{j}); 
        end
        fprintf(stream,'   dxdt = dxdt(:);\n');
        %fprintf(stream,'end\n');
        fprintf(stream,'\n');
    end
    fclose(stream);
    addpath(cd); % add path so that the function is recognized
    if(mode==1 && compile) % compile!
        options = cParserSet('blockSize',1e07);
        ode = [derscript,'.m'];        
        convertToC(indS,ode,{},options);
        compileC([model.name,'_Cvode']);
    end

end