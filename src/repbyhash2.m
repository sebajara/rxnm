function [fstr] = repbyhash2(istr,ca,repl)

% This function is used to take a reaction rate or propensity and
% replace species or parameter symbols by an array construction
%
% [FSTR] = repbyhash(ISTR,CA,REPL)
%
% For example:
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
    fstr = istr;
    sca = sortcabylength(ca,0);
    l = length(ca);
    ind = struct;
    for i = 1:l
        ind.(ca{i}) = i;
    end
    for i = 1:l
        target = sca{i};
        j = ind.(target);
        repby = repl{j};
        fstr = regexprep(fstr,target,repby); % replace strings
    end

end