%% plot spectral peaks and their significance for different upper IPI cutoffs (55ms and 75ms)
cc()
addpath(genpath('src'))
Fs = 10000;
fileNames = {'../spectra/Stern2014_KyriacouManual2017_spec_ipiCutoff75ms_full.mat', '../spectra/Stern2014_FSSStern2014_spec_ipiCutoff75ms_full.mat', '../spectra/Stern2014_FSSCoen2014_spec_ipiCutoff75ms_full.mat';...
   '../spectra/Stern2014_KyriacouManual2017_spec_ipiCutoff55ms_full.mat', '../spectra/Stern2014_FSSStern2014_spec_ipiCutoff55ms_full.mat', '../spectra/Stern2014_FSSCoen2014_spec_ipiCutoff55ms_full.mat'};
label = {'manual Kyriacou et al. (2017)', 'FSS Stern (2014)', 'FSS Coen et al. (2014)'};

% build colormap that is grey for p>0.05
clear cmap
alpha = 0.05;
Cbins = -3:0.1:0; % bins for p-values of peaks
idx = find(10.^Cbins>=0.05, 1, 'first');
cmap(1:idx,:) = limit(summer(idx)-0.05);
cmap(idx+1:length(Cbins),:) = .8;
colormap(cmap);
xlabelPeriod = {'ipi cutoff 75 ms', 'ipi cutoff 55 ms'};

%%
close all
for FREQ = 1:2
   figure('Name', sprintf('Fig. S%d - %s', FREQ+3, xlabelPeriod{FREQ}))
   clf
   for fil = 1:size(fileNames,2) % for each segmentation
      spc = load(fileNames{FREQ, fil});
      token = strsplit(fileNames{FREQ, fil}, {'/' '_'});
      
      dat = load(['../data/' strjoin(token(3:4), '_') '.mat']);
      CSflies = find(contains(dat.flyStrains, 'CS'));
      %% gather all data
      nFlies = length(CSflies);
      Xpeak = nan(nFlies, 1);
      Ppeak = nan(nFlies, 1);
      Ypeak = nan(nFlies, 1);
      Xall = nan(nFlies, 1);
      Pall = nan(nFlies, 1);
      Yall = nan(nFlies, 1);
      for fly= 1:nFlies
         try
            Xpeak(fly) = spc.spec.F{CSflies(fly)}(spc.peak.loc{CSflies(fly)}(argmax(spc.peak.amp{CSflies(fly)}))); % x is save as period in s
            Ppeak(fly) = spc.peak.prob{CSflies(fly)}(argmax(spc.peak.amp{CSflies(fly)}));
            Ypeak(fly) = fly;
            Xall(fly) = spc.spec.F{CSflies(fly)}(spc.peak.loc{CSflies(fly)}(argmax(spc.peak.amp{CSflies(fly)}))); % x is save as period in s
            Pall(fly) = spc.peak.prob{CSflies(fly)}(argmax(spc.peak.amp{CSflies(fly)}));
            Yall(fly) = fly;
         end
      end
      
      [~, Ccnt] = histc(log10(Ppeak), Cbins);
      Ccnt = Ccnt+1;
      subplot(3,1,fil)
      hold on
      for ii = 1:max(Ccnt)
         idx = Ccnt==ii;
         plot(1./(Xpeak(idx)), Ypeak(idx), '.', 'MarkerSize', 22, 'Color', cmap(ii,:));
      end
      if fil==1 % place colorbar in first panel
         colormap(cmap)
         hcb = colorbar('east');
         hcb.Limits = [0 1];
         hcb.Position(3:4) = hcb.Position(3:4)/2;
         title(hcb, 'p-value')
         set(hcb, 'Ticks', linspace(0,1,length(Cbins(1:5:end))), 'TickLabels', round((10.^(Cbins(1:5:end)))*1000)/1000)
      end
      
      if fil==3 % xlabel for last panel only
         xlabel('period [seconds]')
      end
      title(label{fil})
      ylabel('fly#')
      legend('off')
      set(gca, 'XLim', [0 300], 'XTick', 0:50:300)
      set(vline(1./([1/20 1/150])), 'Color', [.6 .6 .6])
      set(vline(1./([1/50 1/70])), 'Color', 'k')
   end
   set(gcas, 'YLim', [0 100])
   
   clp()
end