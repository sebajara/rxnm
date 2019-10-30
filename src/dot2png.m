function [pngfile] = dot2png(dotfile)
    
% This function takes a string representing the file path to a dot
% file (DOTFILE), and uses the dot program from graphviz to create
% the corresponding graph in png format. It return the name of the
% png file. 
% 
% For documentation on dot format see
% http://www.graphviz.org/
%
% [PNGFILE] = dot2png(DOTFILE)
% 
% Sebastian Jaramillo-Riveri
% June, 2013
    
    if(~regmbool(dotfile,'\.dot$'))
        error('Expected filename with .dot extension');
    end
    pngfile = dotfile;
    pngfile = regexprep(pngfile,'\.dot$','.png');
    dotexec = '/usr/bin/dot';
    if(ismac())
        dotexec = '/usr/local/bin/dot';
    end
    [s,r] = system([dotexec,' -Tpng ',dotfile,' -o ',pngfile]);
    if(regmbool(r,'\w'))
        error('Failed making dot graph of the model');
    end
end