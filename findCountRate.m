%bins = parseVideo();

close all
set(groot, 'defaultLegendInterpreter','latex');

counts = bins(:,2);
times = bins(:,1);

exponentialBinsFactor = 1.05;
constantSizedBinsSecondInterval = 3;
logarithmic = 0; % Choose between condensed bins versions


condensedCountsC = zeros(1, floor(times(length(times))/constantSizedBinsSecondInterval));
condensedTimesC = zeros(1, floor(times(length(times))/constantSizedBinsSecondInterval));

lastTimeWritten = -1;
for k = 1:floor(times(length(times))/constantSizedBinsSecondInterval)
    first = find(times > k*constantSizedBinsSecondInterval);
    idx = first(1);
    
    if (times(idx) == lastTimeWritten) % prevent later div-by-0 problems
        continue
    end
    
    condensedCountsC(k) = counts(idx);
    condensedTimesC(k) = times(idx);
    lastTimeWritten = times(idx);
end

condensedCountsC(condensedTimesC == 0) = [];
condensedTimesC(condensedTimesC == 0) = [];

% Create logarithmic intervals.
n = 0;
while (exponentialBinsFactor^n < times(end))
    n = n+1;
end
n = n-1;

condensedCountsL = zeros(1, n);
condensedTimesL = zeros(1, n);
lastTimeWritten = -1;
for k = 1:n
    first = find(times > exponentialBinsFactor^k);
    
    if (first(1) == idx)
        continue
    end
    
    idx = first(1);
    if (times(idx) == lastTimeWritten) % prevent later div-by-0 problems
        continue
    end
    
    condensedCountsL(k) = counts(idx);
    condensedTimesL(k) = times(idx);
    lastTimeWritten = times(idx);
end

condensedCountsL(condensedTimesL == 0) = [];
condensedTimesL(condensedTimesL == 0) = [];

figure
semilogy(condensedTimesC(1:end-1), diff(condensedCountsC)./diff(condensedTimesC),'o');
hold on
set(gca,'TickLabelInterpreter', 'latex');
semilogy(condensedTimesL(1:end-1), diff(condensedCountsL)./diff(condensedTimesL),'x');
xlabel('Time (s)');
ylabel('Count Rate (counts per second)');
t = title('Count Rate vs. Time for Constant and Exponentially Condensed Bins');
set(t, 'Interpreter', 'Latex');
legend('Constant Bins', 'Exponential Bins');

figure
plot(condensedTimesC, condensedCountsC, 'o');
hold on
set(gca,'TickLabelInterpreter', 'latex');
plot(condensedTimesL, condensedCountsL, 'x');
xlabel('Time (s)');
ylabel('Counts');
legend('Constant Bins', 'Exponential Bins', 'Location', 'northwest');
t = title('Counts vs. Time for Constant and Exponentially Condensed Bins');
set(t, 'Interpreter', 'Latex');

% Fit line to late counts.
if (logarithmic)
    condensedTimes = condensedTimesL;
    condensedCounts = condensedCountsL;
    str = 'Exponentially';
    base_str = sprintf('%1.2f', exponentialBinsFactor);
else
    condensedTimes = condensedTimesC;
    condensedCounts = condensedCountsC;
    str = 'Constant';
    base_str = sprintf('%3.0f', constantSizedBinsSecondInterval);
end

first = find(condensedTimes > 180);
cutoff = first(1);

lateTimes = condensedTimes(cutoff:end-1);
lateCountRates = diff(condensedCounts(cutoff:end))./diff(condensedTimes(cutoff:end));

[lateCountRateFit, S1] = polyfit(lateTimes, real(log(lateCountRates)), 1);
lambda = lateCountRateFit(1); % exponential decay constant
A = exp(lateCountRateFit(2))/(-lambda);

figure
errorbar(lateTimes, log(lateCountRates), sqrt(lateCountRates)./lateCountRates, 'x');
hold on
set(gca,'TickLabelInterpreter', 'latex');
plot(lateTimes, polyval(lateCountRateFit, lateTimes, S1));

% uncertainty in parameters of fit:
ste1 = sqrt(diag(inv(S1.R)*inv(S1.R'))./S1.normr.^2./S1.df);

xlabel('Time (s)');
ylabel('Logarithm of Count Rate');
legend('Late Count Rates', 'Linear Fit');
t = title(['Counts vs. Late Time for ' str ' Condensed Bins (Base = ' base_str ')']);
set(t, 'Interpreter', 'Latex');
lambda_str = sprintf('%1.5f', lambda);
lambda_ste_num = ste1(1);
lambda_ste = sprintf('%0.4f', lambda_ste_num);
A_str = sprintf('%4.0f', round(A, -2));
A_ste_num = round(abs(A*(ste1(2) + ste1(1)/lambda)), -2);
A_ste = sprintf('%0.0f', A_ste_num);
exp_str = ['e^{' lambda_str ' \pm ' lambda_ste 't}'];
fit_eq = ['$C_2(t) = ' A_str ' \pm ' A_ste '(1 - ' exp_str ')$'];
x = xlim;
y = ylim;
x_loc = x(1) + .12*(x(2)-x(1));
y_loc = y(1) + .12*(y(2)-y(1));
y_loc2 = y(1) + .07*(y(2)-y(1));
half_life = sprintf('%4.1f', -log(2)/lambda);
hl_ste = sprintf('%4.1f', ste1(1)^2 * log(2)/(lambda*lambda));
t_hl = ['Half life: $' half_life ' \pm ' hl_ste '$ seconds'];
t = text(x_loc, y_loc, fit_eq);
t1 = text(x_loc, y_loc2, t_hl);
set(t, 'Interpreter', 'Latex');
set(t1, 'Interpreter', 'Latex');

% Subtract slow decay counts from early counts.

earlyCutoff = 30;

first = find(condensedTimes > earlyCutoff);
cutoff_idx = first(1);

earlyCounts = condensedCounts(1:cutoff_idx);
condensedTimes;
earlyTimes = condensedTimes(1:cutoff_idx);

subbedEarlyCounts = earlyCounts - A * (1 - exp(lambda * earlyTimes));

earlyCountRates = diff(subbedEarlyCounts)./diff(earlyTimes);

earlyTimes(earlyCountRates < 1) = [];
earlyCountRates(earlyCountRates < 1) = [];

[earlyCountRateFit, S] = polyfit(earlyTimes(1:end-1), real(log(earlyCountRates)), 1);

% uncertainty in parameters of fit:
ste = sqrt(diag(inv(S.R)*inv(S.R'))./S.normr.^2./S.df);

figure
errorbar(earlyTimes(1:end-1), real(log(earlyCountRates)), sqrt(real(log(earlyCountRates)))./real(log(earlyCountRates)), 'x');
hold on
set(gca,'TickLabelInterpreter', 'latex');
plot(earlyTimes, (polyval(earlyCountRateFit, earlyTimes, S)));

lambda2 = earlyCountRateFit(1);
lambda2_ste_num = ste(1);
lambda2_ste = sprintf('%0.4f', lambda2_ste_num);

A2 = exp(earlyCountRateFit(2))/(-lambda2);
xlabel('Time (s)');
ylabel('Logarithm of Count Rate');
legend('Early Count Rates', 'Linear Fit');
t = title(['Counts vs. Early Time for ' str ' Condensed Bins (Base = ' base_str ')']);
set(t, 'Interpreter', 'Latex');
lambda_str = sprintf('%1.5f', lambda2);
A2_str = sprintf('%4.0f', round(A2, -2));
A2_ste_num = round(abs(A2*(ste(2) + ste(1)/lambda2)), -1);
A_ste2 = sprintf('%0.0f', A2_ste_num);

exp_str = ['e^{' lambda_str ' \pm ' lambda2_ste 't}'];
fit_eq = ['$C_1(t) = ' A2_str ' \pm ' A_ste2 '(1 - ' exp_str ')$'];

x = xlim;
y = ylim;
x_loc = x(1) + .15*(x(2)-x(1));
y_loc = y(1) + .15*(y(2)-y(1));
y_loc2 = y(1) + .10*(y(2)-y(1));
half_life = sprintf('%4.1f', -log(2)/lambda2);
hl_ste2 = sprintf('%4.1f', ste(1)^2 * log(2)/(lambda2*lambda2));
t_hl = ['Half life: $' half_life ' \pm ' hl_ste2 '$ seconds'];
t = text(x_loc, y_loc, fit_eq);
t1 = text(x_loc, y_loc2, t_hl);
set(t, 'Interpreter', 'Latex');
set(t1, 'Interpreter', 'Latex');

% Reconstruct plot from two parts

figure
C1 = @(t) A * (1 - exp(lambda * t));
C2 = @(t) A2 * (1 - exp(lambda2 * t));
calcCounts = C1(times) + C2(times);
plot(times, counts);
hold on

plot(times, calcCounts);
plot(times, counts-sqrt(counts), 'c');
C1 = @(t) (A+A_ste_num) * (1 - exp((lambda) * t));
C2 = @(t) (A2+A2_ste_num) * (1 - exp((lambda2) * t));
calcCounts = C1(times) + C2(times);
plot(times, calcCounts, 'r');

plot(times, counts+sqrt(counts), 'c');



C1 = @(t) (A-A_ste_num) * (1 - exp((lambda) * t));
C2 = @(t) (A2-A2_ste_num) * (1 - exp((lambda2) * t));
calcCounts = C1(times) + C2(times);
plot(times, calcCounts, 'r');

set(gca,'TickLabelInterpreter', 'latex');

xlabel('Time (s)');
ylabel('Counts');
legend('Observed Counts', '$C_1(t) + C_2(t)$', 'Uncertainty in Counts', 'Uncertainty in Fit', 'Location', 'northwest');
t = title(['Reconstructed Counts vs. Time for ' str ' Condensed Bins (Base = ' base_str ')']);
set(t, 'Interpreter', 'Latex');
set(gcf, 'Position', [100, 100, 1000, 500]);
