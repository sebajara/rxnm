function [fname] = model2rxnm(model)
% Takes a MODEL struct and creates a file in rxnm format representing
% model symbols, initial conditions, distributions, correlations and
% reactions. It returns the name of the file created.
%
% [FNAME] = model2rxnm(MODEL)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    rxns = fieldnames(model.rxns);
    species = fieldnames(model.species);
    params = fieldnames(model.params);

    nrxns = length(rxns);
    nspec = length(species);
    npara = length(params);

    name = model.name;
    fname = [name,'.rxnm'];
    
    stream = fopen(fname,'w');
    for i = 1:npara
        id = params{i};
        val = model.params.(id).init;
        fprintf(stream,'P: %s %13.5e\n',id,val);
    end
    printLines(stream,{''});

    for i = 1:nspec
        id = species{i};
        val = model.species.(id).init;
        fprintf(stream,'S: %s %13.5e\n',id,val);
    end
    printLines(stream,{''});

    for i = 1:nrxns
        id = rxns{i};
        eq = model.rxns.(id).equation;
        rate = model.rxns.(id).rate;
        fprintf(stream,'R: %s\n   %s\n',eq,rate);
    end
    printLines(stream,{''});

    for i = 1:nspec
        id = species{i};
        if (regmbool(model.species.(id).type,'RV2|RV3'))
            distf = [id,'dist.dat'];
            dist = model.species.(id).dist;
            savematrix(distf,dist);
            fprintf(stream,'D: %s %s\n',id,distf);
        end
    end
    printLines(stream,{''});

    for i = 1:nspec
        id1 = species{i};
        type1 = model.species.(id1).type;
        for j = i:nspec
            if (~(i==j))
                id2 = species{j};
                type2 = model.species.(id2).type;
                if (regmbool(type1,'RV2|RV3') &&...
                    regmbool(type2,'RV2|RV3'))
                    c = model.corr(i,j);
                    fprintf(stream,'C: %s %s %13.5e\n',id1,id2,c);
                end
            end
        end
    end
    printLines(stream,{''});
    fclose(stream);

end