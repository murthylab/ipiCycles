%%
% Fig. 3G
% Influence of IPI threshold on power to detect simulated rhythms
%%

clear all
load('CantonS_KHIPIs_LLR=0.mat')

%%%
%Power analysis
%repeat and calculate % of times get frequency between 50-60 seconds
%(0.0167 - 0.02 Hz)
power = nan(numel(IPI_results),numel(200:50:650));
lomb_results = cell(numel(IPI_results),numel(200:50:650));
thresholded_data = lomb_results;
for sample = 1:numel(IPI_results)
    d = IPI_results(sample).IPI.d;
    t = IPI_results(sample).IPI.t;
    time = 1:t(end);
    fs = 1e4;
    freq = 1/(55*fs);%freq = 1/period
    A = 20; %amplitude 2msec
    x = A *sin(2*pi*freq*t);
    d_sine = x(:) + d(:); %raw data with sine imposed on top
    thresholds = 200:50:650;
    sign = nan(numel(thresholds),1);
    for j = 1:numel(thresholds)
        %keep all IPI values where random # greater than some threshold
        t_thresholded = t(d_sine<thresholds(j));
        d_sine_thresholded = d_sine(d_sine<thresholds(j));
        thresholded_data{sample,j} = [t_thresholded',d_sine_thresholded];
        if numel(d_sine_thresholded) >=1000
            if ~isempty(t_thresholded)%lomb crashes if array empty
                [P,f,alpha] = lomb(d_sine_thresholded,t_thresholded./1e4);
                peak = min(alpha(f>1/60 & f<1/50));
                if peak < 0.05
                    sign(j) = 1;
                else
                    sign(j) = 0;
                end
            end
            lomb_results{sample,j} = [P,f,alpha]; 
        else
            sign(j) = NaN;
            lomb_results{sample,j} = NaN;
        end
    end
    power(sample,:) = sign;
end


%plot Power results
figure(1)
hold on
plot(thresholds/10,nanmean(power),'-o','LineWidth',2)
ylim([0 1.05])
xlim([thresholds(1)/10 thresholds(end)/10])
xlabel('IPI Threshold (msec)','FontSize',26)
ylabel('Power P < 0.05','FontSize',26)
set(gca,'FontSize',24)
set(gca,'box','off')

