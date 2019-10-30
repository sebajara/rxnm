function [status] = printLines(stream,lines)
% Takes a string cell array LINES, and prints to STREAM each value
% as a single line.
%
% [STATUS] = printLines(STREAM,LINES)
%    
% Sebastian Jaramillo-Riveri
% June, 2013    
    
    status = 1;
    for i = 1:length(lines)
        count = fprintf(stream,'%s\n',lines{i});
        if (count == 0) 
            % If it does not print anything, status is turned to 1
            status = 0;
        end
    end

end