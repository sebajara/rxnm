function [fmodel] = cmstoch(imodel)

% Takes a model structure IMODEL, and return a new model FMODEL.
% The new model contains changes related to the stochastic
% interpretation of the model. In particular, the propensety functions.
%
% For now, only changes power functions to a productorial serie. For
% example x^3 -> x*(x-1)*(x-2)
%
% [FMODEL] = cmstoch(IMODEL)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    fmodel = imodel;

    rxns = fieldnames(fmodel.rxns);
    species = fieldnames(fmodel.species);
    params = fieldnames(fmodel.params);

    nrxns = length(rxns);
    nspec = length(species);
    npara = length(params);

    sysep = '(^|\-+|\++|\/+|\(+|\)+|\*+|\^+|$)'; % separator for symbols
                                                 % in a rxn rate function

    [SDM,PDM] = dependmatrix(fmodel);
    SM = stoichmatrix(fmodel);

    % Check for species
    for i = 1:nrxns
        id = rxns{i};
        rate = fmodel.rxns.(id).rate;
        mrate = rate;
        for j = 1:nspec
            if (~(SDM(i,j)==0))
                spc = species{j};
                if(regmbool(mrate,[sysep,spc,'\^\d+']))
                    e = regexprep(mrate,['.*',sysep,spc,'\^(\d+)',sysep],'$2');
                    e = str2num(e);
                    if(e)
                        repby = spc;
                        for q = 2:e
                            repby = [repby,'*(',num2str(1/q),')*(',spc,'-',num2str(q-1),')'];
                        end
                    else
                        error('Error while interpreting a reaction propensity')
                    end
                    repby = ['(',repby,')'];
                    mrate = regexprep(mrate,[spc,'\^\d+'],repby);
                end
            end
        end
        for j = 1:npara
            if (~(PDM(i,j)==0))
                par = params{j};
                if(regmbool(mrate,[sysep,par,'\^\d+']))
                    error('TODO');
                end
            end
        end
        if(~strcmp(mrate,rate))
            fprintf('Reaction propensity of %s was changed:\n',id);
            fprintf('  %s\n',rate);
            fprintf('  ->\n',rate);
            fprintf('  %s\n',mrate);
            fmodel.rxns.(id).rate = mrate;
        end
    end

end