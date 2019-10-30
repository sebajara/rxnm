function [dotfile] = model2dot(model)
% Takes a MODEL structure and creates a digraph in dot format
% representing the relationship between species, parameters, and
% reactions. Stoichiometric relations are represented by full arrows,
% and dependencies on reaction rates are represented by dotted
% arrows. Return the name of the created file.
% 
% [DOTFILE] = model2dot(MODEL)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
    spcc = 'red';
    rxnc = 'blue';
    parc = 'green';
    varc = 'green';
    
    rxns = fieldnames(model.rxns);
    species = fieldnames(model.species);
    params = fieldnames(model.params);
    vars = fieldnames(model.vars);
    
    nrxns = length(rxns);
    nspec = length(species);
    npara = length(params);
    nvars = length(vars);
    
    %arrow = '--';
    %type = 'graph';
    arrow = '->';
    type = 'digraph';
    
    SM = stoichmatrix(model);
    [SDM,PDM,VDM] = dependmatrix(model);
    dotfile = [model.name,'.dot'];
    stream = fopen(dotfile,'w');
    fprintf(stream,'%s %s {\n',type,model.name);
    fprintf(stream,['\n# Species\n']); 
    for i = 1:nspec
        id = species{i};
        name = model.species.(id).name;
        fprintf(stream,' %s [label="%s",color=%s];\n',id,name,spcc);
    end
    fprintf(stream,['\n# Parameters\n']); 
    for i = 1:npara
        id = params{i};
        name = model.params.(id).name;
        fprintf(stream,' %s [label="%s",color=%s,shape=diamond];\n',id,name,parc);
    end
    fprintf(stream,['\n# Variables\n']); 
    for i = 1:nvars
        id = vars{i};
        name = model.vars.(id).name;
        fprintf(stream,' %s [label="%s",color=%s];\n',id,name,varc);
    end
    fprintf(stream,['\n# Reactions\n']); 
    for i = 1:nrxns
        id = rxns{i};
        name = model.rxns.(id).equation;
        fprintf(stream,' %s [label="%s",color=%s,shape=box];\n',id,name,rxnc);
    end
    
    fprintf(stream,['\n# Stoichiometric relations\n']); 
    for i = 1:nrxns
        rxn = rxns{i};
        for j = 1:nspec
            spc = species{j};
            if (SM(i,j)>0) % Is produced by the reaction
                fprintf(stream,' %s %s %s;\n',rxn,arrow,spc);
            elseif(SM(i,j)<0) % Is consumed by the reaction
                fprintf(stream,' %s %s %s;\n',spc,arrow,rxn);
            else
                % is 0, so we do nothing
            end
        end
    end
    fprintf(stream,['\n# Rate dependencies\n']); 
    for i = 1:nrxns
        rxn = rxns{i};
        for j = 1:nspec
            spc = species{j};
            if (~(SDM(i,j)==0)) % Specie is part of the rate function
                fprintf(stream,' %s %s %s [style=dotted];\n',spc,arrow,rxn);
            end
        end
    end
    for i = 1:nrxns
        rxn = rxns{i};
        for j = 1:npara
            par = params{j};
            if (~(PDM(i,j)==0)) % Parameter is part of the rate function
                fprintf(stream,' %s %s %s [style=dotted];\n',par,arrow,rxn);
            end
        end
    end
    for i = 1:nrxns
        rxn = rxns{i};
        for j = 1:nvars
            var = vars{j};
            if (~(VDM(i,j)==0)) % Parameter is part of the rate function
                fprintf(stream,' %s %s %s [style=dotted];\n',var,arrow,rxn);
            end
        end
    end
    fprintf(stream,[' }\n']);
    fclose(stream);
    
end