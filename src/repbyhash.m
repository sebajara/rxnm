function [fstr] = repbyhash(istr,ca,pref,porb,dif)

% This function is used to take a reaction rate or propensity and
% replace species or parameter symbols by an array construction
%
% [FSTR] = repbyhash(ISTR,CA,PREF,PORB,DIF)
%
% For example:
% isrt = 'x*y*z';
% ca = {'x','z','y'};
% repbyhash(isrt,ca,'v',1,0) -> 'v(1)*v(3)*v(2)'
% repbyhash(isrt,ca,'v',0,0) -> 'v[1]*v[3]*v[2]'    
% repbyhash(isrt,ca,'v',0,-1) -> 'v[0]*v[2]*v[1]' 
%    
% TODO: Use varagin for passing PORB and DIF as optionals
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
    fstr = istr;
    sca = sortcabylength(ca);
    l = length(ca);
    ind = struct;
    for i = 1:l
        ind.(ca{i}) = i;
    end
    for i = 1:l
        %target = sca{l-i+1};
        target = sca{i};
        j = ind.(target);
        if (porb)
            repby = [pref,'(',num2str(j+dif),')']; % one has the index
                                                   % bewteen parethesis
        else
            repby = [pref,'[',num2str(j+dif),']']; % and the other
                                                   % bewteen brackets
        end   
        fstr = regexprep(fstr,target,repby); % replace strings
    end

end