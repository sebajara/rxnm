function [fmodel] = rmnonused(imodel)
% It takes a model structure IMODEL and removes unnecessary
% structures.  For now, it only removes reactions that have always
% rates evaluated to 0 (for example if some parameter is 0). It
% decides that by evaluating the rate with all species set to 1, and
% the actual paramater values.
% 
% [FMODEL] = rmnonused(IMODEL)
%
% TODO: check for reachable species and remove non used parameters?
% 
% Sebastian Jaramillo-Riveri
% June, 2013
    
    fmodel = struct;
    fmodel.name = imodel.name;
    fmodel.species = imodel.species;
    fmodel.ind2spc = imodel.ind2spc;
    fmodel.corr = imodel.corr;
    fmodel.params = imodel.params;
    fmodel.ind2par = imodel.ind2par;
    fmodel.sym2ind = imodel.sym2ind;
    fmodel.rxns = struct;

    rxns = fieldnames(imodel.rxns);
    species = fieldnames(imodel.species);
    params = fieldnames(imodel.params);

    inrxns = length(rxns);
    inspec = length(species);
    inpara = length(params);

    SM = stoichmatrix(imodel);
    [SDM,PDM] = dependmatrix(imodel);

    keeprxns = {};
    keepspec = {};
    keeppara = {};
    fnrxns = 0;
    fnspec = 0;
    fnpara = 0;

    % set vector with parameter values
    p = zeros(inpara,1);
    for i = 1:inpara
        par = params{i};
        p(i) = imodel.params.(par).init;
    end

    % decide to keep the rxn based on whether the rate is non zero
    for i = 1:inrxns
        id = rxns{i};
        rxn = imodel.rxns.(id);
        rate = rxn.rate;
        mrate = repbyhash(rate,imodel.ind2par,'p',1,0);
        mrate = repbyhash(mrate,imodel.ind2spc,'x',1,0);
        bool = 0;
        for j = 1:10 % it will try with 10 different sets, if all
                     % evaluate to 10, bool = 0.
            % set vector with species
            x = ones(inspec,1) + j.*rand(inspec,1);
            b = eval(mrate);
            if(b>0)
                bool = 1;
            end
        end
        if(bool)
            fnrxns = fnrxns + 1;
            keeprxns{fnrxns} = id;
        else
            cprintf('err',' Warning:');
            cprintf('text',' 0 rate reaction was removed:');
            cprintf('k',  ' %s\n         %s\n',id,imodel.rxns.(id).equation);
        end
    end
    for i = 1:fnrxns
        id = keeprxns{i};
        fmodel.rxns.(id) = imodel.rxns.(id);
    end
    % 
    
end

