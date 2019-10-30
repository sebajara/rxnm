function [hfname] = model2SSAc(model)
% Takes a MODEL structure and writes down a header file in C,
% implementing a stochastic simulation algorithm for the MODEL.
% Return the file name containing the code.
%    
% Inside that file, you would find two similar functions:
% SSAs and SSAt, with signature
% void SSA (double maxE, double maxt, 
%           double x0[], double p[], 
%           int NTPS, double tps[], 
%           FILE *s){
%  ...
% }
% where maxE is the maximum number of events, maxt is the maximal
% 'bio' time to simulate, x0 are the initial state condition, p are
% the initial parameter condition, NTPS is the number of time points,
% tps are the time points we are interested, and s is a file stream
% where the output will be printed out.
%   
% SSAs will print for each element in tps a new line, containing
% the time, and the state of each specie the same order as they are
% stored. Values are separated by spaces.
%    
% SSAt will print for each element in tps a new line, containing the
% time, the previous time that was printed (or 0 for the first), the
% average value of each specie during that time interval, and then the
% average value of the square of each specie. Values are separated by
% spaces.
%  
% [HFNAME] = model2SSAc(MODEL)
%   
% TODO: I use double sometimes where only int is needed. 
% 
% Sebastian Jaramillo-Riveri
% June, 2013    
    
    rxns = fieldnames(model.rxns);
    species = fieldnames(model.species);
    params = fieldnames(model.params);
    
    nrxns = length(rxns);
    nspec = length(species);
    npara = length(params);
    
    %SM = stoichmatrix(model);
    %[SDM,PDM] = dependmatrix(model);
    
    hfname = [model.name,'_SSA.h'];
    stream = fopen(hfname,'w');
    
    % Definitions for vector constructions
    defl = {['#define NSPE ',num2str(nspec), ' /* Number of species */'],...
            ['#define NPAR ',num2str(npara), ' /* Number of parameters */'],...
            ['#define NRXN ',num2str(nrxns), ' /* Number of reactions */'],...
            ''...
           };
    printLines(stream,defl);
             
    mltadd = {'void v_mltadd(int n, double x1[], double x2[], double s, double x3[]){',...
              ' /* x3 <- x1 + s*x2 */',...
              '  int q = 0;',...
              '  for (q=0;q<n;q++){',...
              '    x3[q] = x1[q] + s*x2[q];',...
              '  }',...
              '}'};
    printLines(stream,mltadd);

    vsum = {'double v_sum (int n, double x[]){',...
            ' /* Sums the elements of x up to n */',...
            '  int q = 0;',...
            '  double sum = 0;',...
            '  for (q=0;q<n;q++){',...
            '    sum += x[q];',...
            '  }',...
            '  return sum;',...
            '}'};
    printLines(stream,vsum);
    
    SSAname = 'SSAs';
    dobewteen = {...
        '    while ((NTPS > ct) && (tps[ct] <= t)){',...
        '      fprintf(s,"%e",tps[ct]);',...
        '      for (q=0;q<NSPE;q++){',...
        '        fprintf(s," %d",(int)x[q]);',...
        '      }',...
        '      fprintf(s,"\n");',...
        '      ct++;',...
        '    }'...
                };
    printSSA (stream,model,SSAname,{''},dobewteen);
    
    % SSA version that prints first and second moments
    SSAname = 'SSAt';
    init = {'  double x1[NSPE]; /* State first moments */',...
            '  double x2[NSPE]; /* State second moments */',...
            '  double xx[NSPE];',...
            '  double it = 0; /* InterIntervals time */'...
           };
    dobewteen = {...
        '    /* Average trajectories over time intervals */',...
        '    if ((NTPS > ct) && (tps[ct] <= t)){',...
        '      /* Add reminder */',...
        '      it += dt-(t-tps[ct]);',...
        '      v_mltadd(NSPE,x1,x,dt-(t-tps[ct]),x1); /* First moment */',...
        '      for (q=0;q<NSPE;q++){',...
        '         xx[q] = x[q]*x[q];',...
        '      }',...
        '      v_mltadd(NSPE,x2,xx,dt,x2); /* Second moment */',....
        '      if (it>0){ /* Normalise by interval length*/',...
        '       for (q=0;q<NSPE;q++){',...
        '         x1[q] = x1[q]/it;',...
        '         x2[q] = x2[q]/it;',...
        '       }',...
        '      }else{',...
        '       for (q=0;q<NSPE;q++){',...
        '         x1[q] = x[q];',...
        '         x2[q] = xx[q];',...
        '       }',...
        '      }',...
        '      fprintf(s,"%e %e",tps[ct]-it,tps[ct]); /* Print initial and final times*/',...
        '      for (q=0;q<NSPE;q++){',...
        '        fprintf(s," %e",x1[q]); /* Print first moments */',...
        '      }',...
        '      for (q=0;q<NSPE;q++){',...
        '        fprintf(s," %e",x2[q]); /* Print second moments */',...
        '      }',...
        '      fprintf(s,"\n");',...
        '      ct++;',...
        '      /* Possibly we cross more than one interval */',...
        '      while((NTPS > ct) && (tps[ct] <= t)){',...
        '       for (q=0;q<NSPE;q++){',...
        '         x1[q] = x[q];',...
        '         x2[q] = xx[q];',...
        '       }',...
        '       fprintf(s,"%e %e",tps[ct-1],tps[ct]); /* Print initial and final times*/',...
        '       for (q=0;q<NSPE;q++){',...
        '         fprintf(s," %e",x1[q]); /* Print first moments */',...
        '       }',...
        '       for (q=0;q<NSPE;q++){',...
        '         fprintf(s," %e",x2[q]); /* Print second moments */',...
        '       }',...
        '       fprintf(s,"\n");',...
        '	ct++;',...
        '      }',...
        '      /* Reset and add reminder */',...
        '       for (q=0;q<NSPE;q++){',...
        '         x1[q] = 0;',...
        '         x2[q] = 0;',...
        '       }',...
        '      it = t-tps[ct-1];',...
        '      v_mltadd(NSPE,x1,x,it,x1); /* First moment */',...
        '      for (q=0;q<NSPE;q++){',...
        '         xx[q] = x[q]*x[q];',...
        '      }',...
        '      v_mltadd(NSPE,x2,xx,it,x2); /* Second moment */',....
        '    }else{',...
        '      it += dt;',...
        '      v_mltadd(NSPE,x1,x,dt,x1); /* First moment */',...
        '      for (q=0;q<NSPE;q++){',...
        '         xx[q] = x[q]*x[q];',...
        '      }',...
        '      v_mltadd(NSPE,x2,xx,dt,x2); /* Second moment */',....
        '    }'...
                };
   
    printSSA (stream,model,SSAname,init,dobewteen);
    fclose(stream);
end

function [] = printSSA (stream,model,SSAname,init,dobewteen)
% The structure of both simulation algorithms are similar, the only
% difference are the name (SSAname), certain lines to instantiate
% extra variables (INIT), and what to do before dobewteen (DOBEWTEEN)
    
    rxns = fieldnames(model.rxns);
    species = fieldnames(model.species);
    params = fieldnames(model.params);
    
    nrxns = length(rxns);
    nspec = length(species);
    npara = length(params);
    
    SM = stoichmatrix(model);
    [SDM,PDM] = dependmatrix(model);
    
    printLines(stream,{['void ',SSAname,' (double maxE, double maxt, double x0[], double p[], int NTPS, double tps[], FILE *s) {']});
    % Variables declarations and initial values
    % '  /* First we make some data type checks */',...
    % '  if (x0 == VNULL || p == VNULL || tps == VNULL) error(E_NULL,"SSAs");',...
    % '  if (x0->dim != NSPE || p->dim != NPAR) error(E_SIZES,"SSAs");',...
    srands = {'  /* Set random seed */',...
              '  struct timeval tv;',...
              '  gettimeofday(&tv, NULL);',...
              '  unsigned char *po = (unsigned char *)&tv;',...
              '  unsigned seed = 0;',...
              '  size_t i;',...
              '  for (i = 0; i < sizeof tv; i++){',...
              '    seed = seed * (UCHAR_MAX + 2U) + po[i];',...
              '  }',...
              '  srand (seed);',...
              ''};
    printLines(stream,srands);
    topl = {...
        '  ',...
        '  /* Variables declarations */', ...
        '  double x[NSPE]; /* Current state */',...
        '  double R[NRXN]; /* Current Propensities */',...
        '  double event = 0; /* Number of events */',...
        '  double t = 0; /* Current time */',...
        '  double rsum = 0; /* For sum of propensities */',...
        '  double lsum = 0; /* For cumsum of propensities */',...
        '  double r1 = 0; /* Random number */',...
        '  double r2 = 0; /* Random number */',...
        '  double dt = 0; /* Time updates */',...
        '  int j, q;',...
        '  int ct = 0;',...
        '  /* Initialise values */',...
        '  for(q=0;q<NSPE;q++){',...
        '    x[q] = x0[q];',...
        '  }'
           };
    printLines(stream,topl);
    
    for i = 1:nrxns
        printPropensities(stream,model,i);
        %fprintf(stream,'  propR%d(R,p,x);\n',i);
    end
    
    % For variables declarations?
    printLines(stream,init);
    
    % Simulation loop
    simtl1 = {...
        '  do {',...
        '    rsum = v_sum(NRXN,R); /* sum of propensities*/',...
        '    r2 = rand() / (double) RAND_MAX; /* rand, 0 and 1*/',...
        '    dt = log(1/r2)/rsum; /* Increase in time */',...
        '    t += dt; /* time update */',...
        };
    printLines(stream,simtl1);
    
    % What to do in between events!!
    printLines(stream,dobewteen);
    
    simtl2 = {...
        '    if (t<maxt){',...
        '      r1 = (rand() / (double) RAND_MAX)*rsum;',...
        '      lsum = 0;',...
        '      for (j=1;j<=NRXN;j++){',...
        '        lsum += R[j-1];',...
        '        if (lsum >= r1){',...
        '          switch(j){'...
             };
    printLines(stream,simtl2);
    
    for i = 1:nrxns
        fprintf(stream,'          case %d: /* %s */\n',i,model.rxns.(rxns{i}).equation);
        printUpState(stream,model,i);
        % fprintf(stream,'            upstateR%d(x);\n',i);
        % update the corresponding propensities
        for q = 1:nrxns
            upd = 0;
            for j = 1:nspec
                if ((~SM(i,j)==0) && (~(SDM(q,j)==0))) 
                    upd = 1;
                end
            end
            if (upd)
                printPropensities(stream,model,q);
                %fprintf(stream,'            propR%d(R,p,x); /* %s */\n',q,model.rxns.(rxns{q}).equation);
            end
        end
        fprintf(stream,'            break;\n');
    end
    simbl = {...
        '          default:',...
        '            break;',...
        '          }',...
        '          break;',...
        '	}',...
        '      }',...
        '      event++;',...
        '    }',...
        '  }while((maxE>event) && (maxt>t));',...
        '  return;',...
        '}',...
        ''...
            };
    printLines(stream,simbl);
end

function [] = printUpState(stream,model,i)
% It prints out the code necessary to update the states given that
% reaction ith was fired.
    rxns = fieldnames(model.rxns);
    species = fieldnames(model.species);
    SM = stoichmatrix(model);
    nrxns = length(rxns);
    nspec = length(species);
    for j = 1:nspec
        spc = species{j};
        ind = model.species.(spc).index;
        c = SM(i,j);
        if (c>0)
            fprintf(stream,'            x[%d] += %d; /* %s */\n',ind-1,c,spc);
        elseif(c<0)
            fprintf(stream,'            x[%d] -= %d; /* %s */\n',ind-1,-1*c,spc);
        end
    end
end

function [] = printPropensities(stream,model,i)
% It prints out the code necessary to update the propensity of
% reaction ith
    rxns = fieldnames(model.rxns);
    nrxns = length(rxns);
    rate = model.rxns.(rxns{i}).rate;
    rate = repbyhash(rate,model.ind2par,'p',0,-1);
    rate = repbyhash(rate,model.ind2spc,'x',0,-1); 
    fprintf(stream,'            R[%d] = %s; /* %s */\n',i-1,rate,model.rxns.(rxns{i}).rate);
end