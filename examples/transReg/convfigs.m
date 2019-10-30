ind = 15;
cudir = cd;
cd 'lowouts/';
load('transReg_low_ss.mat');
x1 = Xss{ind};
cd(cudir);
cd 'highouts1/';
load('transReg_high1_ss.mat');
x2 = Xss{ind};
cd(cudir);
cd 'highouts2/';
model = parserxnm('transReg.rxnm');
simo = makesimo(model,'high2','ss',500000,2287);
[T,X] = loadsimout(simo,{'X'});
x3 = X{1};
cd(cudir);

fig = figure(); hold on;
maxX = max([max(x1) max(x2) max(x3)]);
xbins = [0:maxX];
h1 = myrelfreq(x1,xbins);
h2 = myrelfreq(x2,xbins);
h3 = myrelfreq(x3,xbins);
stairs(xbins,cumsum(h1),'Color','green','LineWidth',2);
stairs(xbins,cumsum(h2),'Color','red','LineWidth',2);
stairs(xbins,cumsum(h3),'Color','blue','LineWidth',2);
set(gca,'XScale','log');
xlabel('X','FontSize',25); ylabel('c.m.f','FontSize',25);
legend({'OFF','ON(1)','ON(2)'});
legend('Location','NorthOutside','Orientation','horizontal');
set(gca,'FontSize',20);
axis([100 20000 0 1]);
imgsize = [0 0 7 7];
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
print(fig,char('cdfs'),'-dpng');

imgsize = [0 0 5 5];
%fig1 = open('highouts1/transReg_high1_ssA_vs_L_vs_X.fig');
%set(fig,'Visible','Off');
%set(fig,'PaperUnits','inches','PaperPosition',imgsize);
%print(fig,char('highouts1/transReg_high2_ssA_vs_L_vs_X.png'),'-dpng');

imgsize = [0 0 5 5];
%fig2 = open('highouts2/transReg_high2_ssA_vs_L_vs_X.fig');
%set(fig,'Visible','Off');
%set(fig,'PaperUnits','inches','PaperPosition',imgsize);
%print(fig,char('highouts2/transReg_high2_ssA_vs_L_vs_X.png'),'-dpng');
    
cudir = cd;
cd 'highouts1/';
model = parserxnm('transReg.rxnm');
simo = makesimo(model,'high1','ss',500000,2287);
[figs] = plotscatter(simo,{{'A','L','X'}});
axis([0 300 0 30 0 1000]);
view(-40,10);
fig1 = figs{1};
cd(cudir);

cd 'highouts2/';
model = parserxnm('transReg.rxnm');
simo = makesimo(model,'high2','ss',500000,2287);
[figs] = plotscatter(simo,{{'A','L','X'}});
axis([0 300 0 30 0 1000]);
view(-40,10);
fig2 = figs{1};
cd(cudir);

imgsize = [0 0 10 10];
set(fig1,'PaperUnits','inches','PaperPosition',imgsize);
print(fig1,char('scatter1.png'),'-dpng');

imgsize = [0 0 10 10];
set(fig2,'PaperUnits','inches','PaperPosition',imgsize);
print(fig2,char('scatter2.png'),'-dpng');