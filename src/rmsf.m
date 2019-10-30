function [sf] = rmsf(si)
% remove space at front and end
    sf = regexprep(si,'^\s*',''); 
    sf = regexprep(sf,'\s*$',''); 
end