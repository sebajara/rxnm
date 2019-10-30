function [model,status] = parserxnm(varargin)
% It takes a FILE with .rxnm extention, parses it, and return a MODEL
% structure. 
%
% Lines begining in # are interpreted as comments.
% 
% Species and Parameters are specified by strings (symbols). They
% need to be unique. Initial values are specified by lines begining
% in 'P:' for parameters, and 'S:' for species. For example
% P: kf 10
% P: kr 1   
% S: A 10
% S: B 10    
% V: l kf*A
% 
% Reactions are specified by a equation (from which the stoichiometry
% is inferred) and a rate function that depends on parameters and
% species. WARNING! Units of parameters and species are not
% checked. It is adviced to think of units of molecules per cell,
% instead of molar; or directly in the discrete setting.
% For example
% R: A + B => AB
%    kf*A*B
% or if reversible
% R: A + B <=> AB
%    kf*A*B
%    kr*AB
% All reversible reactions will be decomposed into two individual
% reactions.
%
% Random variable initial state for species are specified by the
% following kind of line 
% D: A file.dat 
% D: B file.dat 
% where is a text file with 2 columns, the first are species state and
% the second the corresponding mass density.
%
% Correlations are specified as follow
% C: A B 0.8
% which means that the initial state of A and B is correlated with
% value 0.8
% 
% ---------------------------------------------------------------
% 
% The model structure goes as follow
% 
% model.name := string.
% model.species := struct. Contains species structures
% model.params := struct.  Contains parameter structures  
% model.rxns := struct.  Contains reactions structures      
% model.ind2spc := string cell array, contaning species field names
% model.ind2par := string cell array, contaning parameters field names    
% model.ind2var := string cell array, contaning vars field names    
% model.sym2ind := struct. Each parameter and specie field name as
%                  their corresponding index as value.
% 
% Each specie structure goes by
% model.species.(id).name := string (usually the same as field name)
% model.species.(id).index := string (integer > 0)
% model.species.(id).type := string. Either RV1 or RV2. RV1 represents
%                            constant initial condition (see init
%                            attribute), whereas RV2 represents random
%                            variable initial condition (see dist
%                            attribute).
% model.species.(id).init := initial value for RV1 species (also use
%                            for deterministic interpretation, even if
%                            is RV2).
% model.species.(id).dist := n by 2 matrix, where for each row the
%                            first value is a specie state and the
%                            second its mass density.
% 
% Each params structure goes by
% model.params.(id).name := string (usually the same as field name)
% model.params.(id).index := string (integer > 0)
% model.params.(id).type := string. Either RV1 or RV2. RV1 represents
%                           constant initial condition (see init
%                           attribute), whereas RV2 represents random
%                           variable initial condition (see dist
%                           attribute). WARNING! Interpretation of RV2
%                           params are not implemented yet.
% model.params.(id).init := initial value for RV1 species (also use
%                           for deterministic interpretation, even if
%                           is RV2).
% model.params.(id).dist := n by 2 matrix, where for each row the
%                           first value is a specie state and the
%                           second its mass density.
%
% For each reaction, the id is constructed as Rn, where n is a
% counter. The structure goes as follow
% model.rxns.(id).name := string. Not used?
% model.rxns.(id).lhs := string. Left hand side of the reaction equation.
% model.rxns.(id).rhs := string. Right hand side of the reaction equation.
% model.rxns.(id).equation := string. Reaction equation.
% model.rxns.(id).stoich := vector 1 by nspecies, specifying the
%                           stoichiometry of the reaction.
% model.rxns.(id).rate := string. It represents the reaction rate
%                         (or propensety!) written with species and
%                         parameters symbols.
% 
% -------------------------------------------------------------------
%
% [MODEL,STATUS] = parserxnm(FILE)
%
% Sebastian Jaramillo-Riveri
% June, 2013

    endlinecom = '\s*\#.*$';
    
    ext = '\.rxnm$'; % Valid extention for rxnm files
    pline = '^P:';   % regex to identify parameter definitions
    rline = '^R:';   % regex to identify reactions definitions
    sline = '^S:';   % regex to identify species definitions
    dline = '^D:';   % regex to identify distributions definitions
    cline = '^C:';   % regex to identify Correlations definitions
    vline = '^V:';   % regex to identify dynamic variables
    vpline = '^V\+\='; % regex to identify += updates on dynamic variables
    dgline = '^DEFGEN:'; % Generic reactions definition
    gline  = '^GEN:'; % Generic reactions
    
    pvstr = '\s*(\S+)\s+(\S+)\s*.*'; % regex to separate specificiation
                                     % of values
    sdir1 = '^\s*(.*)\s*=>\s*(.*)$'; % regex for unidir reactions
    sdir12 = '^\s*(.*)\s*<=\s*(.*)$'; % regex for unidir reactions
    sdir2 = '^\s*(.*)\s*<=>\s*(.*)$'; % regex for reversible reactions
    sysep = '\-+|\++|\/+|\(+|\)+|\*+|\^+|\s+'; % separator for symbols
                                               % in a rxn rate function
    spsep = '\++|\s+'; % separator for species in a rxn equation

    dval = 0; % default value
    nspec = 0; % current number of species
    npara = 0; % current number of parameters
    nvars = 0; % current number of variables
    nrxns = 0; % current number of reactions
    maxspec = 100; % we need to create a pre made correlation matrix

    % Initiate the structures
    model = struct;
    model.name = '';
    model.species   = struct;
    model.params    = struct;
    model.rxns      = struct;
    model.generic   = struct;
    model.sym2ind   = struct;

    model.ind2spc   = {};
    %model.corr = diag(ones(maxspec,1));
    model.ind2par   = {};
    model.ind2var   = {};

    % Just to make sure is a proper string
    %files = mysplit(ls(varargin),'\s');
    files = varargin;
    
    for f = 1:length(files)
        file = char(files{f});

        if (~ regmbool(file, ext)) % Check extension
            status = 0;
            error(['Wrong file extension: ',file]);
        else
            [fid,mssg] = fopen(file); % Open file
            if (mssg) % Could not open the file...
                error(['File does not exist ',file]);
            else
                model.name = regexprep(file,ext,'');
                lline = fgetl(fid);
                linec = 0;
                lines = {};
                while ischar(lline) % For each line
                    linec = linec+1;
                    lines{linec} = lline;
                    lline = fgetl(fid); % get the next line
                end
                %linec = linec+1;
                %lines{linec} = lline;
                fclose(fid); % Finish parsing the file
                
                parse_lines(lines); % call aux function to parse
                                    % each line
                
            end % finish parsing all files
            rxns = fieldnames(model.rxns); % all rxns keys.
                                           % First we make sure that we have all parameters. Every new 
                                           % symbol from rxn rate will be assumed to be a parameter.
% $$$      % Not in use any more
% $$$             for i = 1:nrxns
% $$$                 rxn = rxns{i};
% $$$                 req = model.rxns.(rxn).rate;
% $$$                 [model,npara] = parserxnreq(model,req,npara);
% $$$             end

            % Now we make some cell arrays to facilitate the mapping
            % bewteen symbols and indexes
            species = fieldnames(model.species); % all species keys
            params = fieldnames(model.params); % all parameters keys
            vars = fieldnames(model.vars); % all variables keys
            model.ind2spc = cell(nspec,1);
            model.ind2par = cell(npara,1); 
            model.ind2var = cell(nvars,1); 
            for i = 1:nspec
                sp = species{i};
                ind = model.species.(sp).index;
                model.ind2spc{ind} = sp;
            end
            for i = 1:npara
                sp = params{i};
                ind = model.params.(sp).index;
                model.ind2par{ind} = sp;
            end
            for i = 1:nvars
                sp = vars{i};
                ind = model.vars.(sp).index;
                model.ind2var{ind} = sp;
            end
            % We will infer stoichiometry from the reactions equations
            for i = 1:nrxns
                rxn = rxns{i};
                eq = model.rxns.(rxn).equation;
                lhs = model.rxns.(rxn).lhs;
                rhs = model.rxns.(rxn).rhs;
                stoich = zeros(1,nspec);
                for j = 1:nspec
                    spc = species{j};
                    ind = model.species.(spc).index;
                    lc = infercoeff(lhs,spc,spsep);
                    lr = infercoeff(rhs,spc,spsep);
                    if (~((lr - lc) == 0))
                        stoich(ind) = lr - lc;
                    end
                    model.rxns.(rxn).stoich = stoich;
                end
            end
            % Correct the correlation matrix
            % model.corr = model.corr(1:nspec,1:nspec);
            if(~validmodel(model))
                error('Invalid model construction');
            end
        end
    end
    
    function parse_lines (llines)
        totallines = length(llines);
        llinecount = 0;
        while(llinecount<totallines)
            llinecount = llinecount+1;
            lline = rmsf(llines{llinecount});
            lline = regexprep(lline,endlinecom,''); % remove comments
            if (regmbool(lline,sline)) % If it specifies a specie
                lline = regexprep(lline,sline,''); % remove sline
                                                   % get symbol and values
                [sym,val] = get_paired_values(lline,pvstr); 
                % Add symbol and value to current species
                [model,nspec] = addsp2struc(model,'species',nspec,sym,val);
            elseif(regmbool(lline,pline)) % If it specifies a parameter
                lline = regexprep(lline,pline,''); % remove pline
                [sym,val] = get_paired_values(lline,pvstr); 
                % Add symbol and value to current parameters
                [model,npara] = addsp2struc(model,'params',npara,sym,val);
            elseif(regmbool(lline,vpline)) % If it specifies a +=variable
                lline = regexprep(lline,vpline,''); % remove vpline
                [sym,val] = get_paired_strings(lline,pvstr);
                % Add symbol and value to current variables
                if(~isfield(model.vars,sym))
                    error(['Incorrect use of +=. Symbol ',sym,'is unknown variable']);
                else
                    model.('vars').(sym).formula = [model.('vars').(sym).formula,'+',val];
                    %[model,nvars] = addsp2struc(model,'vars',nvars,sym,val);
                end
            elseif(regmbool(lline,vline)) % If it specifies a variable
                lline = regexprep(lline,vline,''); % remove vline
                [sym,val] = get_paired_strings(lline,pvstr); 
                % Add symbol and value to current variables
                [model,nvars] = addsp2struc(model,'vars',nvars,sym,val);
            elseif(regmbool(lline,dgline)) % If it defines a generic reaction
                lline = regexprep(lline,dgline,''); % remove gpline
                [sym,val] = get_paired_strings(lline,pvstr); 
                val = regexprep(val,'\(|\)',''); % for now only one
                                                 % parameter is accepted?
                genlines = {};
                linecount = 0;
                while ~(regmbool(lline,'}')) % Add lines until the end
                    linecount = linecount+1;
                    llinecount = llinecount+1;
                    lline = llines{llinecount};
                    %lline = fgetl(fid); % get the next line
                    genlines{linecount} = lline;
                end
                model.generic.(sym) = struct;
                model.generic.(sym).name = sym;
                model.generic.(sym).var  = val;
                model.generic.(sym).code = genlines(1:linecount-1);
            elseif(regmbool(lline,gline)) % If it calls a generic reaction
                lline = regexprep(lline,gline,''); % remove gpline
                [sym,val] = get_paired_strings(lline,pvstr); 
                codelines = model.generic.(sym).code;
                var = model.generic.(sym).var;
                for asdf46 = 1:length(codelines) % replace variable
                                                 % by the argument
                    codelines{asdf46} = regexprep(codelines{asdf46},['_',var],['_',val]);
                end
                parse_lines(codelines);
            elseif(regmbool(lline,rline)) % If it specifies a reaction
                lline = regexprep(lline,rline,''); % remove rline
                rev = 0; % 1 for reversible reactions
                sstr = ''; % string to parse the sides
                if regmbool(lline,sdir2)
                    sstr = sdir2;
                    rev = 1;
                    % get the left hand side
                    lhs = regexprep(lline,sstr,'$1'); 
                    % get the right hand side
                    rhs = regexprep(lline,sstr,'$2');
                elseif regmbool(lline,sdir1)
                    sstr = sdir1;
                    % get the left hand side
                    lhs = regexprep(lline,sstr,'$1'); 
                    % get the right hand side
                    rhs = regexprep(lline,sstr,'$2');
                elseif regmbool(lline,sdir12)
                    sstr = sdir12;
                    % get the left hand side
                    rhs = regexprep(lline,sstr,'$1'); 
                    % get the right hand side
                    lhs = regexprep(lline,sstr,'$2');
                else
                    error('Invalid reaction');
                end
                
                lhs = regexprep(lhs,'^\s+','');
                rhs = regexprep(rhs,'^\s+','');
                lhs = regexprep(lhs,'\s+$','');
                rhs = regexprep(rhs,'\s+$','');
                % reform the equation in one direction
                eq = [lhs,' => ',rhs];
                % increase number of reactions
                nrxns = nrxns + 1;
                % make a dummy id
                id = ['R',num2str(nrxns)];
                % add new species from the sides
% $$$                 % Not in use any more
% $$$                 [model,nspec] = parserxnside(model,lhs,nspec);
% $$$                 [model,nspec] = parserxnside(model,rhs,nspec);
                model.rxns.(id) = struct; % initialize struct
                model.rxns.(id).name = ''; % name?
                model.rxns.(id).equation = eq; %equation
                model.rxns.(id).lhs = lhs; % 
                model.rxns.(id).rhs = rhs; % 
                llinecount = llinecount+1;
                lline = llines{llinecount};
                %%lline = fgetl(fid); % next line shoud have the equation
                model.rxns.(id).rate = rmsf(lline); % put the rate equation
                if (rev) % redo for the other direccion if reversible
                    eq = [rhs,' => ',lhs];
                    nrxns = nrxns + 1;
                    id = ['R',num2str(nrxns)];
                    model.rxns.(id) = struct;
                    model.rxns.(id).name = '';
                    model.rxns.(id).equation = eq;
                    model.rxns.(id).lhs = rhs; % 
                    model.rxns.(id).rhs = lhs; % 
                    llinecount = llinecount+1;
                    lline = llines{llinecount};
                    %%lline = fgetl(fid); % next line shoud have the equation
                    model.rxns.(id).rate = rmsf(lline);
                end
            elseif(regmbool(lline,dline)) % If it specifies a
                                          % distribution
                lline = regexprep(lline,dline,''); % remove dline
                [sym,val] = get_paired_values(lline,pvstr);
                if (isfield(model.species,sym))
                    model.species.(sym).type = 'RV2';
                    model.species.(sym).dist = load_table(val);
                end
            elseif(regmbool(lline,cline)) % If it specifies a
                                          % correlation
                lline = regexprep(lline,cline,''); % remove cline
                error('Correlations not supported yet\n');
            end
        end
    end
    
% $$$     function [fmodel,nf] = parserxnside(imodel,side,ni)
% $$$     % Go troughs all symbols in a reaction equation side and create
% $$$     % structures for new species.
% $$$         fmodel = imodel;
% $$$         nf = ni;
% $$$         if (length(side))
% $$$             csps = mysplit(side,spsep); % split to find species
% $$$             for q = 1:length(csps)
% $$$                 csp = csps{q};
% $$$                 if(length(csp) && (~regmbool(csp,'^\d+$')))
% $$$                     % if we have a value and is not a number
% $$$                     if(isfield(model.params,csp)||isfield(model.vars,csp))
% $$$                         error(['In a reaction equation, Param/Var "',csp,'" is being used as a specie']);
% $$$                     end
% $$$                     if(~isfield(model.species,csp))
% $$$                         error(['In a reaction equation, Symbol "',csp,'" is unknown specie']);
% $$$                     end
% $$$                     %[fmodel,nf] = addsp2struc(fmodel,'species',nf,csp,dval);
% $$$                 end
% $$$             end
% $$$         end
% $$$     end
% $$$     
% $$$     function [fmodel,nf] = parserxnreq(imodel,req,ni)
% $$$     % this function reads a reaction equation, and find all symbols
% $$$     % Then creates params structures for each new symbol.
% $$$         fmodel = imodel;
% $$$         nf = ni;
% $$$         if(length(req))
% $$$             csps = mysplit(req,sysep); %FIX THIS?
% $$$             for q = 1:length(csps)
% $$$                 csp = csps{q};
% $$$                 if(length(csp) ...
% $$$                    && (~regmbool(csp,'^\d+$'))... % not a number
% $$$                    && (~isfield(model.sym2ind,csp)))
% $$$                     error(['In: ',eq,'. Unknown symbol "',csp,'".']);
% $$$                     %[fmodel,nf] = addsp2struc(fmodel,'params',nf,csp,dval);
% $$$                 end
% $$$             end
% $$$         end
% $$$     end




end