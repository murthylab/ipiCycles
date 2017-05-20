%Fig.2
%2A Left Panel - Original Data
clear all
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

%2A Middle Panel: 10 min sine
clf
time = 1:t(end);
fs = 1e4;
f = 1/(55*fs);%freq = 1/period
A = 20; %amplitude 2msec
x = A *sin(2*pi*f*t);

plot(t/60/1e4,x/10,'LineWidth',6)
xlim([0 10])
ylim([-30 30])
set(gca,'box','off')
set(gca,'FontSize',24)
xlabel('Time (min)','FontSize',26)
title('IPI cycle','FontSize',30)


%2A Right Panel: raw data with sine imposed on top
clf
d_sine = x(:) + d(:); 
yy = smooth(t/60/1e4,d_sine,200,'rloess');
plot(t/60/1e4,d_sine/10,'.k')
hold on
plot(t/60/1e4,yy/10,'Color', [0.617 0.14 0.56],'LineWidth',6);
ylim([10 70])
xlim([0 10])
set(gca,'box','off')
set(gca,'FontSize',24)
xlabel('Time (min)','FontSize',26)
title('Original Data + IPI Cycle','FontSize',30)

%Fig.2B: original data
clf
yy = smooth(t/60/1e4,d,200,'rloess');
plot(t/60/1e4,d/10,'.k')
ylim([20 65])
xlim([0 45])
hold on
set(gca,'box','off')
ylabel('Inter-Pulse Interval (ms)','FontSize',26)
set(gca,'FontSize',24)
xlabel('Time (min)','FontSize',26)
plot(t/60/1e4,yy/10,'Color', [0.617 0.14 0.56],'LineWidth',3);
title('Original Data','FontSize',30)
hold off

%Fig.2C: LS of original data
clf
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

%Fig.1D: Original data + IPI cycle
clf
yy_sine = smooth(t/60/1e4,d_sine,200,'rloess');
plot(t/60/1e4,d_sine/10,'.k')
ylim([20 65])
xlim([0 45])
hold on
set(gca,'box','off')
ylabel('Inter-Pulse Interval (ms)','FontSize',26)
set(gca,'FontSize',24)
xlabel('Time (min)','FontSize',26)
plot(t/60/1e4,yy_sine/10,'Color', [0.617 0.14 0.56],'LineWidth',3);
title('Original Data + IPI Cycle','FontSize',30)
hold off

%Fig.1E: LS of Original data + IPI cycle
clf
[P,f,alpha] = lomb(d_sine,t/1e4);
[a,z] = significance(d_sine,t);
figure(3)
plot(1./f,P,'k')
xlim([20 150])
styles = {':','-.','--'};
%a = [0.001 0.01 0.05];
for i = 1:numel(styles)
    line([20,130],[z(i),z(i)],'Color','k','LineStyle',styles{i});
    text(131,z(i),strcat('\alpha = ',num2str(a(i))),'fontsize',16); 
end
xlabel('Period (sec)','FontSize',26)
ylabel('Power','FontSize',26)
set(gca,'FontSize',24)
set(gca,'box','off')

% Fig.2F: 5% data (675 IPIs) + IPI cycle 
% Note that data are selected randomly each time and plot will 
% therefore look different each time code is run
clf
j = 0.05;
num = numel(d_sine);
rnd_num = rand(num,1);
t_thresholded = t(rnd_num<j);
d_sine_thresholded = d_sine(rnd_num<j);
yy_one = smooth(t_thresholded/60/1e4,d_sine_thresholded,15,'rloess');
plot(t_thresholded/60/1e4,d_sine_thresholded/10,'.k')
ylim([20 65])
xlim([0 45])
hold on
set(gca,'box','off')
ylabel('Inter-Pulse Interval (ms)','FontSize',26)
set(gca,'FontSize',24)
xlabel('Time (min)','FontSize',26)
plot(t_thresholded/60/1e4,yy_one/10,'Color', [0.617 0.14 0.56],'LineWidth',3);
title('5% Original Data','FontSize',30)
hold off

% Fig.2G: LS of 5% data (675 IPIs) + IPI cycle 
% Note that data are selected randomly each time and plot will 
% therefore look different each time code is run

clf
[P,f,alpha] = lomb(d_sine_thresholded,t_thresholded/1e4);
[a,z] = significance(d_sine_thresholded,t_thresholded);
figure(3)
plot(1./f,P,'k')
xlim([20 150])
styles = {':','-.','--'};
%a = [0.001 0.01 0.05];
for i = 1:numel(styles)
    line([20,130],[z(i),z(i)],'Color','k','LineStyle',styles{i});
    text(131,z(i),strcat('\alpha = ',num2str(a(i))),'fontsize',16); 
end
xlabel('Period (sec)','FontSize',26)
ylabel('Power','FontSize',26)
set(gca,'FontSize',24)
set(gca,'box','off')
