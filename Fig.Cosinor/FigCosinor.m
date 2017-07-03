cc()
d{1} = load('Stern2014_KyriacouManual2017_spec_ipiCutoff75ms_original');
d{2} = load('Stern2014_KyriacouManual2017_spec_ipiCutoff75ms_shuffled');
% d{3} = load('../spectra/Stern2014_KyriacouManual2017_spec_ipiCutoff55ms_full.mat');
%%
clf
for ii = 1:length(d)
   mySubPlot(length(d)+1,2,ii,1)
   plot(1./d{ii}.a.F, log2(1./d{ii}.a.spec))
   set(gca, 'XLim', [20 150], 'YLim', [10 18])
end

for ii = 1:length(d)
   mySubPlot(length(d)+1,2,ii,2)
   hold on
   
   for fly = 1:length(d{ii}.peak.significant)
      try
         plot(1./d{ii}.peak.significant{fly},fly, '.')
      end
   end
   set(gca, 'XLim', [20 150], 'YLim',[0 length(d{ii}.peak.significant)])
end
%%
clear cnt;
bins = 20:5:150;
for ii = 1:length(d)
   cnt(:,ii) = hist(1./vertcat(d{ii}.peak.significant{:}), bins);
end
mySubPlot(length(d)+1,2,length(d)+1,1)
plot(bins, log2(cnt), 'LineWidth', 2)
ylabel('log_2 count')
mySubPlot(length(d)+1,2,length(d)+1,2)
plot(bins, cnt, 'LineWidth', 2)
ylabel('linear count')
legend({'original', 'shuffled'})

   