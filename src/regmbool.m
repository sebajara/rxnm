function [boolean] = regmbool(str,match)

% This function takes a string STR and a regular expression MATCH and
% return 1 or 0 depending on whether the regular expression match the
% string or not.
%
% [BOOLEAN] = regmbool(STR,MATCH)
% 
% Sebastian Jaramillo-Riveri
% June, 2013
    
    if (char((regexp(str, match, 'match'))))
        boolean = 1;
    else
        boolean = 0;
    end

end