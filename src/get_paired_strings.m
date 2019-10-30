function [s1,s2] = get_paired_strings(l,pvstr)
% This function returns two strings as specified by the regex string
% pvstr. Also, it makes sure that the first and second value are not
% equal.
    s1 = regexprep(l,pvstr,'$1'); % get the symbol
    s2 = regexprep(l,pvstr,'$2'); % get the value
    if (strcmp(s1,s2))
        error('Error while using get_paired_strings');
    end

end