function [cam,pos] = cellregexm(ca,regex)

% Takes an string cell array CA and a regular expression string
% REGEX. Return all the elements that match the REGEX (CAM) and the
% positions of those elements POS.
%
% [CAM,POS] = cellregexm(CA,REGEX)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    pos = 0;
    nmatchs = 0;
    for i = 1:length(ca)
        if (regmbool(ca{i},regex))
            nmatchs = nmatchs + 1;
            pos(nmatchs) = i;
        end
    end
    if(pos)
        cam = ca{pos};
    else
        cam = {''};
    end
end