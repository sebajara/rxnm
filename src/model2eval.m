function [evalscript] = model2eval(model)
% Takes a MODEL struct and creates a Matlab 
%
% [EVALSCRIPT] = model2eval(MODEL)
%
% Sebastian Jaramillo-Riveri
% Aug 2014
    
    rxns = fieldnames(model.rxns);
    species = fieldnames(model.species);
    nrxns = length(rxns);
    nspec = length(species);
    npara = length(fieldnames(model.params));
    nvars = length(fieldnames(model.vars));
    SM = stoichmatrix(model);
    
    evalscript = [model.name,'_eval'];
    
    stream = fopen([evalscript,'.m'],'w');
    
    fprintf(stream,'function [result] = evalinmodel (p,ts,xs,expr)\n');
    fprintf(stream,'%%  %s\n',model.name);

    fprintf(stream,'\n');
    for i = 1:npara
        fprintf(stream,'   %s = p(%s);\n',model.ind2par{i},num2str(i));
    end
    fprintf(stream,'\n');
    fprintf(stream,'[npoints,~] = size(xs);\n');
    fprintf(stream,'result = zeros(npoints,1);\n');
    fprintf(stream,'\n');
    fprintf(stream,'for asdf = 1:npoints\n');
    fprintf(stream,'  x = xs(asdf,:);\n');
    fprintf(stream,'  t = ts(asdf);\n');
    for i = 1:nspec
        fprintf(stream,'     %s = x(%s);\n',model.ind2spc{i},num2str(i));
    end
    fprintf(stream,'\n');
    for i = 1:nvars
        formula = model.vars.(model.ind2var{i}).formula;
        fprintf(stream, '    %s = %s;\n',model.ind2var{i},formula);
    end
    fprintf(stream,'\n');
    
    for j = 1:nspec
        fprintf(stream, '     dx_%s = 0',species{j});
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
    
    fprintf(stream,'  result(asdf) = eval(expr);\n');
    fprintf(stream,'  end\n');
    fprintf(stream,'end\n');
    fprintf(stream,'\n');
    fclose(stream);
end