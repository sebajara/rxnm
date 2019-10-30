function [status] = printInclude (stream,libs)
% Prints each of the elements of LIBS to STREAM formated as C
% libraries (include)
%
% [STATUS] = printInclude(STREAM,LIBS)
% 
% Sebastian Jaramillo-Riveri
% June, 2013    
status = 1;
for i = 1:(length(libs))
    count = fprintf(stream,'#include %s\n',libs{i});
    if (count == 0) 
        % If it does not print anything, status is turned to 1
        status = 0;
    end
end

end