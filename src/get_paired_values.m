function [s,v] = get_paired_values(l,pvstr)
% This function returns two values as specified by the regex string
% pvstr. Also, it makes sure that the first and second value are not
% equal.
    s = regexprep(l,pvstr,'$1'); % get the symbol
    v = regexprep(l,pvstr,'$2'); % get the value
    if (strcmp(s,v))
        error('Error while using get_paired_values');
    elseif(regmbool(v,'^\d'))
        v = str2num(v); % return a proper number
    else
        % ??
    end

end
