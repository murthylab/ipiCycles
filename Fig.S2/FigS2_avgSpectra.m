% plot average IPI spectra for the different annotations
addpath('../code')
cc()
fileNames = {'Stern2014_KyriacouManual2017_spec_ipiCutoff75ms_full.mat', ...
             'Stern2014_FSSStern2014_spec_ipiCutoff75ms_cut.mat',...
             'Stern2014_FSSCoen2014_spec_ipiCutoff75ms_cut.mat'};
%%
clf
for fil = 1:length(fileNames)
   load(['../spectra/' fileNames{fil}], 'a');
   a.flyNames(cellfun(@isempty, a.flyNames)) = {''};
   CSidx = contains(a.flyNames, 'CS');
   if fil==1
      a.flyNames(CSidx)
   end
   sum(CSidx)
   subplot(3,1,fil)
   h1 = plot(1./a.F, a.spec(:,CSidx), 'Color', [.6 .6 .6]);
   hold on
   h2 = plot(1./a.F, nanmean(a.spec(:,CSidx),2), 'k', 'LineWidth', 1.5);
   if fil == 1
      legend([h1(1) h2], {'individual fly', 'average'}, 'Box', 'off')
   end
   axis('tight')
   ylabel('Power')
end
set(gcas, 'XScale', 'log', 'XTick', [20 50:50:150], 'YLim', [0 15])
xlabel('Period (s)')
clp()
