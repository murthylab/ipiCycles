% compares the manually segmented data from Kyriacou et al. to results obtained by the FlySongSegmenter with 
% different parameters: the conservative ones used by Stern et al. (2014)
% and those used by Coen et al. (2014)
% plot detection rates as a function of 1. pulse amplitude and 2. pulse
% carrier frequency as well as histograms for the number of consecutively missed pulses
addpath('../code')
cc()
segmentNames = {'Kyriacou', 'David','Pip'};
%% load pulse shapes and 
p = load('../data/Stern2014_KyriacouManual2017_pulses.mat','recId', 'pulseTimes', 'pulseShapes');
p.pulseShapes = double(p.pulseShapes)/1000;     % pulse shapes are save as uint16 after scaling up 1000x to save space - we reverse that here
pulseCarrier = getPulseFreq(p.pulseShapes');    % get pulse carrier frequency...
pulseAmplitude = max(abs(p.pulseShapes),[],2);  % ...and pulse amplitude
%% load pulse times for the three segmenters (Kyriacou manual, David FSS, Pip FSS)
dat{1} = load('../data/Stern2014_KyriacouManual2017.mat');  % manual pulse times from Kyriacou et al. (2017)
dat{2} = load('../data/Stern2014_FSSStern2014.mat');           % automatically segmented with parameters from Stern et al. (2014)
dat{3} = load('../data/Stern2014_FSSCoen2014.mat');             % or Coen et al. (2014)
%% compare pulse times from the three different segmentations for all recordings
tolerance = 5/1000; % jitter (in ms) within which pulse are scored as the same
allPulseAmps = cell(2,1);
allPulseCars = cell(2,1);
for fly = 1:length(dat{1}.flyNames) % for each recording
   %%
   disp(dat{1}.flyNames{fly})
   % pulse times, amplitudes and carrier frequencies
   tmp = dat{1}.pulseTimes(:,fly);
   pulseTimes{fly, 1} = tmp(~isnan(tmp));
   pulseAmps{fly} = pulseAmplitude(p.recId==fly,1);
   pulseCars{fly} = pulseCarrier(p.recId==fly,1);
   disp([size(pulseTimes{fly, 1},1), size(pulseAmps{fly},1)])
%    % keep only unique pulses - deal with pulses that were called more than once during manual segmentation
   [~, uniIdx] =  uniquetol(pulseTimes{fly,1}, tolerance/1000*2);
   pulseTimes{fly,1} = pulseTimes{fly,1}(uniIdx);
   pulseAmps{fly} = pulseAmps{fly}(uniIdx);
   pulseCars{fly} = pulseCars{fly}(uniIdx);
   disp([size(pulseTimes{fly, 1},1), size(pulseAmps{fly},1)])
   idx(1,fly) = fly;
   %%
   for dt = 1:length(dat)-1
      % find current Kyriacou fly in the other data sets
      idx(dt+1,fly) = find(strcmp(dat{dt+1}.flyNames, dat{1}.flyNames{fly}));
      % cut pulseTimes to Kyriacou times
      tmp = dat{dt+1}.pulseTimes(:,idx(dt+1,fly));
      tmp(tmp<min(pulseTimes{fly, 1})-tolerance | tmp>max(pulseTimes{fly, 1})+tolerance) = [];
      pulseTimes{fly, dt+1} = tmp(~isnan(tmp));
   end
   % get errors for all pulses
   [confMat{fly}, eventMat, pulseId, pulseT, pulseGroup] = idPulses(pulseTimes(fly,:), tolerance);
  
   %%
   detectVec = eventMat(logical(eventMat(:,1)),2); % only keep true positives
   missOnset = find(diff(detectVec)==-1);
   for ms = 1:length(missOnset)-1
      consecutivePulsesMissed{fly}(ms) = find(detectVec(missOnset(ms)+1:end)==1,1,'first');
   end
   %%
   % get amp of missed/called pulses
   for ii = 1:2
      pulseAmpsCalled(:, ii, fly) = grpstats(pulseAmps{fly}(1:sum(eventMat(:,1)==1)), eventMat(eventMat(:,1)==1,ii+1), 'mean');
      pulseCarsCalled(:, ii, fly) = grpstats(pulseCars{fly}(1:sum(eventMat(:,1)==1)), eventMat(eventMat(:,1)==1,ii+1), 'mean');
      allPulseAmps{ii} = [allPulseAmps{ii}; [pulseAmps{fly}(1:sum(eventMat(:,1)==1)), eventMat(eventMat(:,1)==1,ii+1)]];
      allPulseCars{ii} = [allPulseCars{ii}; [pulseCars{fly}(1:sum(eventMat(:,1)==1)), eventMat(eventMat(:,1)==1,ii+1)]];
   end
end
%% merge and normalize confusion matrices to get FP/TP rates
cmAll = zeros(size(confMat{1}));
for fly = 1:length(dat{1}.flyNames)
   cmAll = cmAll + confMat{fly};
end
cmaTP = squeeze(cmAll(2,:,1,2:3))'./squeeze(sum(cmAll(2,:,1,2:3),2));      % true positives - all pulses
cmaFP = squeeze(cmAll(1,2,1,2:3))./sum(squeeze(cmAll(:,2,1,2:3)))';        % false positives - all pulses
%% segmenter performance (TP rate) as function pulse amplitude and carrier frequency
cols = lines(2);
figure('Name', '4 - performance of FlySongSegmenter')
for seg = 1:2
   subplot(1,6,1)
   hold on
   bar(seg, cmaTP(seg,2), 0.5, 'EdgeColor', 'none', 'FaceColor', cols(seg,:))
   ylabel('true positive rate')
   
   subplot(1,6,2)
   hold on
   bar(seg, cmaFP(seg,1), 0.5, 'EdgeColor', 'none', 'FaceColor', cols(seg,:))
   ylabel('false positive rate')
end
subplot(1,6,1)
set(gcas, 'XTick', 1:2, 'XTickLabel', {'Stern (2014)', 'Coen et al. (2014)'},...
   'YLim', [0 1], 'YTick', 0:0.2:1, 'XLim', [0.2 2.8])

mySubPlot(2,3,1,2)
hold on
for seg = 1:2
   [xx, yy, xxE, yyE] = tfAdapt(allPulseAmps{seg}(:,1), allPulseAmps{seg}(:,2), 64);
   plot(xx,yy)
end
ylabel('detection probability')
set(gcls, 'LineWidth', 2)
set(gca, 'XLim', [0 1.15], 'YLim', [0 1])
legend({'Stern (2014)', 'Coen et al. (2014)'},  'Box','off', 'Location','SouthEast')

mySubPlot(2,3,2,2)
[cnt, bin] = hist(allPulseAmps{seg}(:,1), linspace(0, 1.2, 64));
plot(bin, normalizeSum(cnt), 'k', 'LineWidth', 2)
xlabel('pulse amplitude [V]')
ylabel('probability')
set(gca, 'XLim', [0 1.15], 'YLim', [0 0.15])

% freq of missed pulses
mySubPlot(2,3,1,3)
hold on
for seg = 1:2
   [xx, yy, xxE, yyE] = tfAdapt(allPulseCars{seg}(:,1), allPulseCars{seg}(:,2), 64);
   plot(xx,yy);
   vline(300)
end
set(gcls, 'LineWidth', 2)
ylabel('detection probability')
set(gca, 'XLim', [0 700], 'YLim', [0 1])

mySubPlot(2,3,2,3)
cla
hold on
[cnt0, bin] = hist(allPulseCars{seg}(:,1), linspace(0, 1000, 64));
plot(bin, normalizeSum(cnt0), 'k', 'LineWidth', 2)
xlabel('pulse carrier frequency [Hz]')
ylabel('number of pulses')
set(gca, 'XLim', [0 700], 'YLim', [0 0.15])

clp()

seg=1;% manual segmenter
% soft pulses missed out of all pulses missed
softPulses = allPulseAmps{seg}(:,1)<0.05;
disp(sum(allPulseAmps{seg}(softPulses,2)==0)/sum(allPulseAmps{seg}(:,2)==0))
% "fast" pulses missed out of all pulses missed
fastPulses = allPulseCars{seg}(:,1)>300;
disp(sum(allPulseCars{seg}(fastPulses,2)==0)/sum(allPulseCars{seg}(:,2)==0))

%% plot number of consecutively missed pulses
figure('Name', '3A - consecutively missed IPIs')
stretchMissed = [consecutivePulsesMissed{:}]-1;
[cnt,bin] = hist(stretchMissed, 1:40);
stretchMissedNumberOfPulses = bin.*cnt; % stretch of given length times number of occurrence of stretches of that length
CDF = cumsum(stretchMissedNumberOfPulses)./sum(stretchMissedNumberOfPulses);

bar(bin,CDF, 0.8, 'EdgeColor', 'none', 'FaceColor', [.6 .6 .6])
hold on
bar(1,CDF(1), 0.8, 'EdgeColor', 'none', 'FaceColor', [1 .3 .3])
vline(1.5)
title(sprintf('%d percent of all missed pulses are one-skips', round(100*CDF(1))))
xlabel('skip size')
ylabel('CDF of missed pulses in skip events of size N')
set(gca, 'YLim', [0 1])
set(gcas, 'XLim', [min(bin)-1 max(bin)+1], 'XTick', [1 5:5:max(bin)])

clp()
