cc()
randLabel = {'original', 'shuffled'};
for ii = 1:length(randLabel)
   d{ii} = load(['Stern2014_KyriacouManual2017_spec_ipiCutoff75ms_' randLabel{ii}]);
end
%%
cols = lines(2);
clf
for ii = 1:length(d)
   mySubPlot(length(d)+1,2,ii,1)
   plot(1./d{ii}.a.F, log2(1./d{ii}.a.spec), 'Color', [.6 .6 .6])
   hold on
   plot(1./d{ii}.a.F, log2(mean(1./d{ii}.a.spec,2)), 'Color', cols(ii,:), 'LineWidth', 1.5)
   set(gca, 'XLim', [20 150], 'YLim', [10 18])
   title(randLabel{ii})
end
ylabel('log_2 cosinor amplitude')
xlabel('period [s]')

for ii = 1:length(d)
   mySubPlot(length(d)+1,2,ii,2)
   hold on
   for fly = 1:length(d{ii}.peak.significant)
      try
         plot(1./d{ii}.peak.significant{fly},fly, '.', 'Color', cols(ii,:))
      end
   end
   set(gca, 'XLim', [20 150], 'YLim',[0 length(d{ii}.peak.significant)])
   title(randLabel{ii})
end
xlabel('fly #')
xlabel('period [s]')

%%
clear cnt;
bins = 10:10:160;
for ii = 1:length(d)
   cnt(:,ii) = hist(1./vertcat(d{ii}.peak.significant{:}), bins);
end
mySubPlot(length(d)+1,2,length(d)+1,1)
plot(bins, cnt, 'LineWidth', 2)
ylabel('linear count')
xlabel('period [s]')
legend(randLabel)
mySubPlot(length(d)+1,2,length(d)+1,2)
plot(bins, log2(cnt), 'LineWidth', 2)
ylabel('log_2 count')
xlabel('period [s]')
set(gcas, 'XLim', [20 150])

clp()
if mfilename()
   figexp(mfilename(), 1.0, 1.0)   
end