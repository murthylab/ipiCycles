%%
% Fig 3b
% Histograms of all IPI values from automatically annotated and manually annotated 
% data from Canton-S recordings
%%

load('/Users/sternd/Dropbox (HHMI)/Writing/KH cycle rebuttal/fromPaperAnalyzedData/>1000IPIs/CantonS_KHIPIs_LLR=0.mat')
all_ipi = [];
for i = 1:numel(IPI_results)
    all_ipi= cat(2,all_ipi,IPI_results(i).IPI.d);
end

load('/Users/sternd/Dropbox (HHMI)/Writing/KH cycle rebuttal/Kyriacou_manual_data/KyriacouManual_data_Stern.mat')
edges = [15:1:75];
histogram(ipi(:,15:39).*1e3,edges,'Normalization','probability','FaceColor','none')
xlim([15 75])
hold on
histogram(all_ipi./10,edges,'Normalization','probability')
hold off