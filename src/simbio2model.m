function [model] = simbio2model(simbiomodel)
% This function will take a simbiology model structure (SIMBIOMODEL),
% and will return a model in our internal format. 
%
% [MODEL] = simbio2model(SIMBIOMODEL)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
    rmvchars = '\s+|\-';
    
    % Initiate the structures
    model = struct;
    model.name = regexprep(get(simbiomodel,'Name'),rmvchars,'_');
    model.species = struct;
    model.params = struct;
    model.rxns = struct;
    model.sym2ind = struct;
    
    nspec = length(simbiomodel.Species);
    model.ind2spc = cell(nspec,1);
    for i = 1:nspec
        spc = simbiomodel.Species(i);
        id = regexprep(get(spc,'Name'),rmvchars,'_');
        if(~isfield(model.species,id))
        model.species.(id).name = get(spc,'Name');
        model.species.(id).index = i;
        model.species.(id).init = get(spc,'InitialAmount');
        model.species.(id).type = 'RV1';
        model.ind2spc{i} = id;
        else
            cprintf('err',' Warning\n');
        end
    end
    model.corr = diag(ones(nspec,1));
    
    npara = 0;
    for i = 1:length(simbiomodel.Parameters)
        par = simbiomodel.Parameters(i);
        npara = npara + 1;
        id = regexprep(get(par,'Name'),rmvchars,'_');
        if(~isfield(model.params,id))
            model.params.(id).name = get(par,'Name');
            model.params.(id).index = npara;
            model.params.(id).init = get(par,'Value');
            model.params.(id).type = 'RV1';
        else
            cprintf('err',' Warning\n');
        end
    end
    
    nrxns = 0;
    for i = 1:length(simbiomodel.Reactions);
        nrxns = nrxns + 1;
        rxn = simbiomodel.Reactions(i);
        id = ['R',num2str(nrxns)];
        model.rxns.(id).name = get(rxn,'Name');
        eq = get(rxn,'Reaction');
        eq = regexprep(eq,'null','');
        rev = get(rxn,'Reversible');
        if (rev)
            eq = regexprep(eq,'<->','=>');
        else
            eq = regexprep(eq,'->','=>');
        end
        rate = get(rxn,'ReactionRate');
        kl = get(rxn,'KineticLaw');
        params = get(kl,'Parameters');
        for j = 1:(length(params))
            par = params(j);
            pid = regexprep(get(par,'Name'),rmvchars,'_');
            if(~isfield(model.params,pid))
                npara = npara + 1;
                model.params.(pid).name = get(par,'Name');
                model.params.(pid).index = npara;
                model.params.(pid).init = get(par,'Value');
                model.params.(pid).type = 'RV1';
                model.ind2par{npara} = pid;
            elseif(model.params.(pid).init == get(par,'Value'))
                cprintf('err',' Warning:');
                cprintf('text',' Attempt to redifine parameter');
                cprintf('k',' %s',pid);
                cprintf('text',' in reaction ');
                cprintf('-comment','%s (%s)\n',id,model.rxns.(id).name);
                cprintf('text','          no action taken, it has the same value.\n');
            else
                pid2 = [pid,'_',id];
                cprintf('err',' Warning:');
                cprintf('text',' In reaction ');
                cprintf('-comment','%s (%s)\n',id,model.rxns.(id).name);
                cprintf('text','          parameter ');
                cprintf('k','%s',pid);
                cprintf('text',' was modified to ');
                cprintf('k','%s',pid2);
                cprintf('text','.\n');
                npara = npara + 1;
                model.params.(pid2).name = get(par,'Name');
                model.params.(pid2).index = npara;
                model.params.(pid2).init = get(par,'Value');
                model.params.(pid2).type = 'RV1';
                model.ind2par{npara} = pid2;
                rate = regexprep(rate,pid,pid2);
            end
        end
        stoich = zeros(1,nspec);
        simstoich = get(rxn,'Stoichiometry');
        prods = get(rxn,'Products');
        react = get(rxn,'Reactants');
        for j = 1:length(react)
            sid = regexprep(react(j).Name,rmvchars,'_');
            if(regmbool(react(j).Name,rmvchars))
                eq = regexprep(eq,['\[',react(j).Name,'\]'],sid);
                rate = regexprep(rate,['\[',react(j).Name,'\]'],sid);
            end
            ind = model.species.(sid).index;
            stoich(ind) = simstoich(j);
        end
        for j = 1:length(prods)
            sid = regexprep(prods(j).Name,rmvchars,'_');
            if(regmbool(prods(j).Name,rmvchars))
                eq = regexprep(eq,['\[',prods(j).Name,'\]'],sid);
                rate = regexprep(rate,['\[',prods(j).Name,'\]'],sid);
            end
            ind = model.species.(sid).index;
            stoich(ind) = simstoich(j+length(react));
        end
        model.rxns.(id).equation = eq;
        if (rev)
            rates = mysplit(rate,'-');
            model.rxns.(id).rate = rates{1};
            if(~regmbool(rate,'-') || (length(rates)>2))
                fprintf('%s\n',rate);
                error('Could not infer how to separate forward and reverse rate functions');
            end
        else
            model.rxns.(id).rate = rate;
        end
        model.rxns.(id).stoich = stoich;
        if (rev)
            model.rxns.(id).name = [model.rxns.(id).name,' fow'];
            nrxns = nrxns + 1;
            id = ['R',num2str(nrxns)];
            model.rxns.(id).name = [get(rxn,'Name'), ' rev'];
            model.rxns.(id).equation = regexprep(eq,'=>','<=');
            model.rxns.(id).rate = rates{2};
            model.rxns.(id).stoich = -1*stoich;
        end
    end
    params = fieldnames(model.params);
    model.ind2par = cell(npara,1);
    for i = 1:npara
        id = params{i};
        ind = model.params.(id).index;
        model.ind2par{ind} = id;
    end
    
    % we are going to save it in rxnm and reparse just in case
    [mfile] = model2rxnm(model);
    model = parserxnm(mfile);
    model2rxnm(model);
end