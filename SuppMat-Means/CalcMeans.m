% Calculate mean and std for automated and manually annotated IPIs

%per0
clear all
load('per0_KHIPIs_LLR=0.mat')
all_ipi = [];
for i = 1:numel(IPI_results)
    all_ipi= cat(2,all_ipi,IPI_results(i).IPI.d);
end
mean(all_ipi)/10
% 41.0214
std(all_ipi)/10
% 8.0877

%perL
clear all
load('perL_KHIPIs_LLR=0.mat')
all_ipi = [];
for i = 1:numel(IPI_results)
    all_ipi= cat(2,all_ipi,IPI_results(i).IPI.d);
end
mean(all_ipi)/10
% 37.6274
std(all_ipi)/10
% 5.9926

%perS
clear all
load('perS_KHIPIs_LLR=0.mat')
all_ipi = [];
for i = 1:numel(IPI_results)
    all_ipi= cat(2,all_ipi,IPI_results(i).IPI.d);
end
mean(all_ipi)/10
% 40.8109
std(all_ipi)/10
% 6.7182

% simNJ 
clear all
load('simNJ_KHIPIs_LLR=0.mat')
all_ipi = [];
for i = 1:numel(IPI_results)
    all_ipi= cat(2,all_ipi,IPI_results(i).IPI.d);
end
mean(all_ipi)/10
% 43.1422
std(all_ipi)/10
% 9.1256

% CantonS 
clear all
load('CantonS_KHIPIs_LLR=0.mat')
all_ipi = [];
for i = 1:numel(IPI_results)
    all_ipi= cat(2,all_ipi,IPI_results(i).IPI.d);
end
mean(all_ipi)/10
% 34.3526
std(all_ipi)/10
% 7.4576

% Manually annotated perL
clear all
load('KyriacouManual_data.mat')
all_ipi = [];
for i = 1:14
    all_ipi= cat(1,all_ipi,ipi(:,i));
end
nanmean(all_ipi)*1000
% 37.5050
nanstd(all_ipi)*1000
% 6.5880

% Manually annotated CantonS
clear all
load('KyriacouManual_data.mat')
all_ipi = [];
for i = 15:39
    all_ipi= cat(1,all_ipi,ipi(:,i));
end
nanmean(all_ipi)*1000
% 33.4322
nanstd(all_ipi)*1000
% 7.1103



