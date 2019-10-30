function [] = savematrix(file,M)

% Saves a matrix (M) into a file (FILE) by printing each row into a
% line, and each column value separated by a space. Prints everything
% with float precision.
%
% savematrix(FILE,M)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
    stream = fopen(file,'w');

    [n,m] = size(M);

    for i = 1:n
        fprintf(stream,'%f',M(i,1));
        if(m>1)
            for j = 2:m
                fprintf(stream,' %f',M(i,j));
            end
        end
        fprintf(stream,'\n',M(i,1));
    end
    fclose(stream);
    
end