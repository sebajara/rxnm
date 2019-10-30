function [fname] = procsim(simo)

% This function creates the C code necessary to do simulation
% 'SIM' structure is as follow:
%    
%   simo.model := a reaction model structure
%   simo.ts    := a vector. time points.
%   simo.type  := string. 'st' for trajectories, and 'ss' for
%                states. 
%   simo.nrep  := integer. Number of repetitions.
%   simo.name  := string. Mnemonic name.
%   
% This function will generate a file by the name of 
%   <simo.model.name>_<simo.name>_<simo.type>.c
% which will contain the implementation of the simulation.  Also will
% create the SSA implementation of the model by calling model2SSAc.
%
% TYPES: 
% ss :=> 
% The printed C code would create single output file containing the
% state of species at each of the specified times during simulation,
% one line per state.
% st :=> 
% The printed C code would create one output file per simulation. On
% it, it would print the average of each state, and the average of the
% square of each state during the time intervals specified by TS.
% 
% [FNAME] = procsimo(SIMO)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
name = simo.name;
model = simo.model;
type = simo.type;
ts = simo.ts;
nrep = simo.nrep;

rxns = fieldnames(model.rxns);
species = fieldnames(model.species);
params = fieldnames(model.params);

nrxns = length(rxns);
nspec = length(species);
npara = length(params);
ntps = length(ts);

maxE = 10000000000; % For now default, doubt we will need something higher
maxt = max(ts);

model2SSAc(model);

% type of job cheeking
if (strcmp(type,'ss'))
    fname = [model.name,'_',name,'_ss.c'];
elseif(strcmp(type,'st'))
    fname = [model.name,'_',name,'_st.c'];
else
    error('Invalid simulation type for procsimo');
end
stream = fopen(fname,'w');
% Libraries and some functions definitions
libs = {'<stdio.h>','<stdlib.h>','<math.h>',...
        '<time.h>','<limits.h>',...
        ['"',model.name,'_SSA.h"']...
       };
printInclude(stream,libs);
printLines(stream,{''});

% Definitions for vector constructions
defl = {['#define NSPE ',num2str(nspec), ' /* Number of species */'],...
        ['#define NPAR ',num2str(npara), ' /* Number of parameters */'],...
        ['#define NRXN ',num2str(nrxns), ' /* Number of reactions */'],...
        ['#define NTPS ',num2str(ntps), ' /* Number of time points */'],...
        ['#define NREP ',num2str(nrep), ' /* Number of replicas */'],...
        ''...
       };
printLines(stream,defl);

topl = {'int main () {',...
        '  double x0[NSPE]; /* Initial State */',...
        '  double xt0[NSPE]; /* Temporal initial state */',...
        '  double p[NPAR]; /* Parameters */',...
        '  double tps[NTPS]; /* Time points */',...
        '',...
        ['  double maxE = ',num2str(maxE),'; /* Maximal number of events */'],...
        ['  double maxt = ',num2str(maxt),'; /* Maximal time */'],...
        '  double rd = 0;',...
        '  int i = 0, qd = 0;',...
        '  FILE *ofp;  char outputF[255];',...
        '',...
       };
printLines(stream,topl);
printLines(stream,{'  /* Time points */'});
for i = 1:length(ts)
    fprintf(stream,'  tps[%d] = %13.5e;\n',i-1,ts(i));
end

% declare new extra variables, such as distributions
for i = 1:nspec
    id = species{i};
    na = model.species.(id).name;
    typ = model.species.(id).type;
    ind = model.species.(id).index;
    if(strcmp(typ,'RV2') || strcmp(typ,'RV3'))
        dist = model.species.(id).dist;
        cdist = cumsum(dist(:,2));
        [nv,~] = size(dist);
        fprintf(stream,'  double dX%d[%d][2];/* c.m.f for initial state of %s */\n',ind,nv,na);
        for j = 1:nv
            fprintf(stream,'  dX%d[%d][0] = %13.5e;dX%d[%d][1] = %13.5e;\n',ind,j-1,dist(j,1),ind,j-1,cdist(j));
        end
        % we could convert distributions into a vector of the inverse of the
        % cumulative distribution function, over equally distanced
        % over the interval [0,1] (fminuspdf)
        % evals{evc} = fminuspdf(model.species.(id).dist);
        % 
        % for j = 1:length(evals{evc})
        %    fprintf(stream,'  dX%d[%d] = %13.5e;\n',ind,j-1,evals{evc}(j));
        % end
        % esistrs = [esistrs,',dX',num2str(ind)];
    end
end


% set of all constant initial species and parameters (RV1)
printinit(stream,model,'RV1');
srands = {'  /* Set random seed */',...
          '  struct timeval tv;',...
          '  gettimeofday(&tv, NULL);',...
          '  unsigned char *po = (unsigned char *)&tv;',...
          '  unsigned seed = 0;',...
          '  size_t z;',...
          '  for (z = 0; z < sizeof tv; z++){',...
          '    seed = seed * (UCHAR_MAX + 2U) + po[z];',...
          '  }',...
          '  srand (seed);',...
          ''};
printLines(stream,srands);

% If we print states (type == 'ss') we want one output file
if (strcmp(type,'ss'))
    % output file for the simulations
    fprintf(stream,'  sprintf(outputF,"%s_%s_ss.out");\n',model.name,name);
    printLines(stream,{'  ofp = fopen(outputF, "w");'});
end

%% simulation loop!!!!
printLines(stream,{'  for(i=0;i<NREP;i++){'});
% during each iteration we will set the initial conditions that are
% random variables (RV2|RV3)
printinit(stream,model,'RV2|RV3');
if (strcmp(type,'ss')) 
    flines = {'    SSAs(maxE,maxt,x0,p,NTPS,tps,ofp);',...
              '   }',...
              '  fclose(ofp);',...
              '}'...
             };
    printLines(stream,flines);
elseif (strcmp(type,'st'))
    % If we print trajectories (type == 'st') we want one output file per simulation
    fprintf(stream,'    sprintf(outputF,"%s_%s_st%%d.out",i);\n',model.name,name);
    flines = {'    ofp = fopen(outputF, "w");',...
              '    SSAt(maxE,maxt,x0,p,NTPS,tps,ofp);',...
              '    fclose(ofp);',...
              '   }',...
              '}'...
             };
    printLines(stream,flines);
    fclose(stream);
end

end

function [] = printinit(stream,model,vartype)
% this function prints the code to set the initial values
    species = fieldnames(model.species);
    params = fieldnames(model.params);
    nspec = length(species);
    npara = length(params);
    
    xinit = 'x0';
    pinit = 'p';
    
    % initial species
    for i = 1:nspec
        id = species{i};
        name = model.species.(id).name;
        typ = model.species.(id).type;
        ind = model.species.(id).index;
        if (strcmp(typ,'RV1') && regmbool(typ,vartype))
            val = model.species.(id).init;
            fprintf(stream,'  %s[%d] = %13.5e; /* %s */\n',xinit,ind-1,val,name);
        elseif((strcmp(typ,'RV2') || strcmp(typ,'RV3')) ...
               && regmbool(typ,vartype))
            dist = model.species.(id).dist;
            fprintf(stream,'    rd = rand() / (double) RAND_MAX;\n');
            [nv,~] = size(dist);
            fprintf(stream,'    for(qd=0;qd<%d;qd++){\n',nv);
            fprintf(stream,'       if(dX%d[qd][1] >= rd){\n',ind);
            fprintf(stream,'         xt0[%d] = dX%d[qd][0]; /* %s */\n',ind-1,ind,name);
            fprintf(stream,'         break;\n');
            fprintf(stream,'       }\n');
            fprintf(stream,'    }\n');
        end
    end
    
    if(regmbool('RV2|RV3',vartype))
        C = model.corr;
        
        % Check if we have correlations bewteen different initial
        % species states
        cc = zeros(1,nspec);
        nc = 1;
        for i = 1:nspec
            bool = 0;
            for j = 1:nspec
                if(~(i==j) && (C(i,j)>0))
                    bool = 1;
                end
            end
            if(bool)
                cc(i) = nc;
                nc = nc + 1;
            end
        end
        if (sum(cc)) % correct in case of correlations
            ccp = find(cc>0);
            [xs,stds] = getspcmean(model);
            Q = C(ccp,ccp).*[stds(ccp)'*stds(ccp)];
            [L,p] = chol(Q,'lower'); % Cholesky factorization.
            if (p == 0) % successful factorization
                for i = 1:nspec
                    id = species{i};
                    ind = model.species.(id).index;
                    typ = model.species.(id).type;
                    if ((strcmp(typ,'RV2') || strcmp(typ,'RV3')) && ...
                        regmbool(typ,vartype) && ...
                        cc(ind))
                        % to correct for correlations we do
                        % <x_i> + floor(L*((x - <x>)./std(x)))
                        fprintf(stream,'    %s[%d] = %d+floor(0',xinit,ind-1,round(xs(ind)));
                        for j = 1:nspec
                            id2 = species{j};
                            ind2 = model.species.(id2).index;
                            if(cc(ind2))
                                if(stds(ind2)>0)
                                    l = L(cc(ind),cc(ind2));
                                    if(l>0)
                                        fprintf(stream,'+%13.5e*((xt0[%d]-%13.5e)/%13.5e)',l,ind2-1,xs(ind2),stds(ind2));
                                    end
                                    else
                                    error(['Zero standard deviation for initial specie distribution']);
                                end
                            end
                        end
                        fprintf(stream,');\n');
                    end
                end
            else
                cprintf('err',' Warning:');
                cprintf('text',' Unsuccessful Cholesky Factorization. Correlations not included for simulations.\n');...
            end
        else
            for i = 1:nspec
                    id = species{i};
                    ind = model.species.(id).index;
                    typ = model.species.(id).type;
                if ((strcmp(typ,'RV2') || strcmp(typ,'RV3')) && ...
                        regmbool(typ,vartype))
                fprintf(stream,'    %s[%d] = xt0[%d];\n',xinit,ind-1,ind-1);   
                end
            end
        end
    end
    % set parameter values
    for i = 1:npara
        id = params{i};
        name = model.params.(id).name;
        typ = model.params.(id).type;
        ind = model.params.(id).index;
        if (strcmp(typ,'RV1') && regmbool(typ,vartype))
            val = model.params.(id).init;
            fprintf(stream,'  %s[%d] = %13.5e; /* %s */\n',pinit,ind-1,val,name);
        elseif(strcmp(typ,'RV2') || strcmp(typ,'RV3'))
            error('Random Variables for Parameters is not implemented yet');
        end
    end
    
end