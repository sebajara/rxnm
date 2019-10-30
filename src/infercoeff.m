function [coeff] = infercoeff(side,spc,spsep)
% Get the coefficient accompanying a specie in the side of a rxn
% equation. 
    coeff = 0;
    if (regmbool(side,['(^|',spsep,')',spc,'(',spsep,'|$)']))
        mstr = ['(\d+?)\s*','(?:',spc,')'];
        c = regexp(side,mstr,'tokens');
        [~,nc] = size(c);
        if (~isempty(c))
            for w = 1:nc
                coeff = coeff+str2num(char(c{w}));
            end
        else
            coeff = 1;
        end
    end
end
