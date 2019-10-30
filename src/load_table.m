function [data] = load_table (filename,varargin)

% Parse a delimited text file FILENAME containing numerical values.
% 
% [DATA] = load_table(FILENAME, 'OptionString', optionvalue ... )
% 
% Options: 
% 'Separator' := string. Regular expression delimiting
%                values. Default '\s+|,'.
% 'Comment' := string. Regular expression for discarting lines
%              Default: '^#'.
% 'NewLine' := string. Regular expression for the end of the line
%              Default: '\n'.
% 
% Sebastian Jaramillo-Riveri 
% June, 2013
   
    sep = '(?:\s+|,)';
    comm = '^#';
    newl = '\n';
    skipt1st = 0;
    
    nop = length(varargin);
    if (nop>0)
        if(mod(nop,2)==0)
            for i = 1:(nop/2)
                op = varargin{i}; % option string
                val = varargin{i+1}; % value
                if (strcmp(op,'Separator'))
                    sep = val;
                elseif(strcmp(op,'Comment'))
                    comm = val;
                elseif(strcmp(op,'NewLine'))
                    newl = val;
                elseif(strcmp(op,'Skipt1st'))
                    skipt1st = val;
                else
                    error('Invalid option string');
                end
            end
        else
            error('Options should be paired as option-string and value');
        end
    end
    
    if (~exist(char(filename), 'file'))
        error('File does not exists');
    end
    
    fid=fopen(char(filename));
    % the end of line might be system dependent
    lines = textscan(fid,'%s','HeaderLines',0,'Delimiter',newl); 
    fclose(fid);
    lines = lines{1};

    [nlines,~] = size(lines);
    % we should infer the number of columns by the first line
    ncols = length(mysplit(lines{1},sep));

    data = zeros(nlines-1,ncols);
    nvals = 0;
    
    st = 1;
    if(skipt1st)
        st = 2;
    end
    
    for i = st:nlines
        l = lines{i};
        if ((~strcmp(l,'')) || (~regmbool(l,comm)))
            nvals = nvals + 1;
            vals = mysplit(l,sep);
            if(length(vals)<ncols)
                error(['not enough columns at line %d',i])
            end
            for j = 1:ncols
                data(nvals,j) = str2num(vals{j});
            end
        end
    end
    % resize the output in case of comment lines
    data = data(1:nvals,1:ncols);
end