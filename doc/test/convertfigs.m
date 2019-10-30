imgsize = [0 0 5 8];
fig = open('outs1/rna1_doc1_det.fig');
set(fig,'Visible','Off');
subplot(3,1,1);
axis([0 100000 0 11]);
subplot(3,1,2);
axis([0 100000 0 5]);
subplot(3,1,3);
axis([0 100000 0 175]);
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
print(fig,char('outs1/rna1_doc1_det.png'),'-dpng');
close(fig);

fig = open('outs1/rna1_doc1_st.fig');
set(fig,'Visible','Off');
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
print(fig,char('outs1/rna1_doc1_st.png'),'-dpng');
close(fig);

imgsize = [0 0 12 10];
fig = open('outs1/rna1_doc1_ss.fig');
set(fig,'Visible','Off');
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
print(fig,char('outs1/rna1_doc1_ss.png'),'-dpng');
close(fig);

imgsize = [0 0 5 5];
fig = open('outs2/rna2_doc2_ssG10_vs_G20.fig');
set(fig,'Visible','Off');
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
axis(gca,[0 21 0 31]);
print(fig,char('outs2/rna2_doc2_ssG10_vs_G20.png'),'-dpng');
close(fig);

fig = open('outs3/rna2_doc3_ssG10_vs_G20.fig');
set(fig,'Visible','Off');
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
print(fig,char('outs3/rna2_doc3_ssG10_vs_G20.png'),'-dpng');
close(fig);

fig = open('outs2/rna2_doc2_ssM1_vs_M2.fig');
set(fig,'Visible','Off');
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
axis(gca,[0 280 0 510]);
print(fig,char('outs2/rna2_doc2_ssM1_vs_M2.png'),'-dpng');
close(fig);

fig = open('outs3/rna2_doc3_ssM1_vs_M2.fig');
set(fig,'Visible','Off');
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
print(fig,char('outs3/rna2_doc3_ssM1_vs_M2.png'),'-dpng');
close(fig);

imgsize = [0 0 10 10];
fig = open('outs4/rna2_doc4_scanss_2.fig');
set(fig,'Visible','Off');
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
print(fig,char('outs4/rna2_doc4_scanss_2.png'),'-dpng');
close(fig);

fig = open('outs4/rna2_doc4_scandet_2.fig');
set(fig,'Visible','Off');
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
print(fig,char('outs4/rna2_doc4_scandet_2.png'),'-dpng');
close(fig);

imgsize = [0 0 7 12];
fig = open('outs5/rna1_doc5_scandet_2.fig');
set(fig,'Visible','Off');
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
print(fig,char('outs5/rna1_doc5_scandet_2.png'),'-dpng');
close(fig);

imgsize = [0 0 20 20];
fig = open('outs5/rna1_doc5_scanss_2.fig');
set(fig,'Visible','Off');
set(fig,'PaperUnits','inches','PaperPosition',imgsize);
for i = 1:4
    subplot(2,2,i);
    set(gca,'FontSize',25);
    xlabel('k01','FontSize',30); ylabel('k10','FontSize',30);
    if(i==1)
        zlabel('<G1>','FontSize',30);
    elseif(i==2)
        zlabel('CV(G1)','FontSize',30);
    elseif(i==3)
        zlabel('<M>','FontSize',30);
    elseif(i==4)
        zlabel('CV(M)','FontSize',30);
    end
end
print(fig,char('outs5/rna1_doc5_scanss_2.png'),'-dpng');
close(fig);