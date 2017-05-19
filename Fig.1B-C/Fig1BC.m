%Fig 1B
load('CantonS_KHIPIs_LLR=0.mat')
d = IPI_results(14).IPI.d;
t = IPI_results(14).IPI.t;
yy = smooth(t/60/1e4,d,200,'rloess');
plot(t/60/1e4,d/10,'.k')
ylim([10 70])
xlim([0 10])
hold on
set(gca,'box','off')
ylabel('Inter-Pulse Interval (ms)','FontSize',26)
set(gca,'FontSize',24)
xlabel('Time (min)','FontSize',26)
plot(t/60/1e4,yy/10,'Color', [0.617 0.14 0.56],'LineWidth',6)
hold off

%Fig 1C
clf
d = d(t<10*60*1e4);%first ten minutes of data
t = t(t<10*60*1e4);
[P,f,alpha] = lomb(d,t/1e4);
[a,z] = significance(d,t/1e4);
plot(1./f,P,'k')
hold on
xlim([20 150])
styles = {':','-.','--'};
for i = 1:numel(styles)
    line([20,130],[z(i),z(i)],'Color','k','LineStyle',styles{i});
    text(131,z(i),strcat('\alpha = ',num2str(a(i))),'fontsize',16); 
end
xlabel('Period (sec)','FontSize',26)
ylabel('Power','FontSize',26)
set(gca,'FontSize',24)
set(gca,'box','off')
hold off
