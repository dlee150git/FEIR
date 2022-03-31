close all; clear; clc

% read data line by line as string array
filename = 'sweep-4496.4nm-58mW-830nm-100mW-PMT10-2.5MHz-1MHz-loop256-new-sample-reverse.txt';
data = readmatrix(filename);
wavenumbers = readmatrix([filename(1:end-4), '_P.txt']);

% find number of spectra in data
nData = length(data);
nSweeps = 0;
t = 1000;

for j = 1:nData
    t1 = data(j,1);
    if t1 < t 
        nSweeps = nSweeps+1;
    end
    t = t1;
end

data_full = permute(reshape(data,nData/nSweeps,[],3),[1,3,2]);

% Input sample parameters
dye = 'ATTO647N 100uM DMSO';
sample_prep = '10mm &25mm CaF2 sandwich';
objective = '1.05NA 2mm obj/Ge lens';
bond = 'triple';
param = '(';
pks1 = zeros(nSweeps,1);
bg1 = zeros(nSweeps,1);
noise1 = zeros(nSweeps,1);
pks2 = zeros(nSweeps,1);
bg2 = zeros(nSweeps,1);
noise2 = zeros(nSweeps,1);
x_master = zeros(nSweeps,1);
power = zeros(nSweeps,1);

for i =1:nSweeps
    data_slice = data_full(:,:,i);
    x = data_slice(:,1); y1 = -data_slice(:,2); y2 = data_slice(:,3);
    [~, low] = max(y2);
    [~, high] = max(y2);
    low = low-10;
    high = high+10;

    f1 = polyfit(x([1:low, high:end]), y1([1:low, high:end]), 1);
    f2 = polyfit(x([1:low, high:end]), y2([1:low, high:end]), 1);
    trend1 = polyval(f1, x);
    trend2 = polyval(f2, x);
    y1Out = y1 - trend1;
    y2Out = y2 - trend2;
    xq = linspace(x(1),x(end),length(x));
    y1Out = interp1(x, y1Out, xq, 'spline');
    y2Out = interp1(x, y2Out, xq, 'spline');
    %% Plot raw data and baseline-corrected data
    figure
    subplot(2,2,1)
    hold on
    plot(x, y1, '.-','Linewidth',1)
    plot(x, trend1, ':', 'Linewidth',1)
    legend('DC Raw data', 'Trend')

    subplot(2,2,2)
    hold on
    plot(x, y2, '.-','Linewidth',1)
    plot(x, trend2, ':', 'Linewidth',1)
    legend('AC Raw data', 'Trend')
    
    %% Process data to get peak heights, S/N and S/B ratios
    subplot(2,2,3)
    threshold1 = max(y1Out)/1.25;
    findpeaks(y1Out(2:end), x(2:end), 'MinPeakHeight', threshold1, 'Annotate','extents','WidthReference','halfheight')
    xlabel('Stage position (mm)'); ylabel('PMT readings (V)')
    [pks_temp, locs, w] = findpeaks(y1Out, 'MinPeakHeight', threshold1);
    if ~isempty(pks_temp)
        pks1(i) = max(pks_temp);
    else
        pks1(i) = max(y1Out);
    end
    bg1(i) = mean(trend1);
    noise1(i) = std(y1Out([1:low, high:end]));
        
    subplot(2,2,4)
    threshold2 = max(y2Out)/1.25;
    findpeaks(y2Out(2:end), x(2:end), 'MinPeakHeight', threshold2, 'Annotate','extents','WidthReference','halfheight')
    xlabel('Stage position (mm)'); ylabel('AC signal (V)')
    [pks_temp, locs, w] = findpeaks(y2Out, 'MinPeakHeight', threshold2);
    if ~isempty(pks_temp)
        pks2(i) = max(pks_temp);
    else
        pks2(i) = max(y2Out);
    end
    
    bg2(i) = mean(trend2);
    noise2(i) = std(y2Out([1:low, high:end]));
    
    sample_str = sprintf('%s-%s-%s',dye,sample_prep, objective);
% %     filename_str = strrep(str,'_','-');
%     if strcmp(bond,'triple')
%         IV_str = sprintf('%s %s (%.2f cm-1)',iv_str,IV,1/str2double(iv_str)*1e7);
%     else 
%         IV_str = sprintf('%s %s (%.2f cm-1)',iv_str,IV,(2/str2double(iv_str)-1/1031.2)*1e7);
%     end
    IV_str = sprintf('%.2f cm-1', wavenumbers(i,1));
%     
    result_str = sprintf('AC Peak height=%.3f bg=%.3f sb=%.3f',pks2(i),bg2(i),pks2(i)/bg2(i));
    
    sgtitle([IV_str result_str])

end
%%
sb1 = pks1./bg1;
sn1 = pks1./noise1;
sb2 = pks2./bg2;
sn2 = pks2./noise2;

figure(1000)
subplot(2,2,1)
% plot(x_master, pks2.*max(power)./power, '.-')
% plot(x_master, pks2./power, '.-')
plot(wavenumbers(:,1), pks2, '.-')
title('Peak height')
subplot(2,2,2)
% plot(x_master, sb2./power, 'o-')
plot(wavenumbers(:,1), sb2, 'o-')
title('S/B ratio')
subplot(2,2,3)
% plot(x_master, sn2./power, 'o-')
plot(wavenumbers(:,1), sn2, 'o-')
title('S/N ratio'); xlabel('Wavenumbers cm-1');
subplot(2,2,4)
% plot(x_master, bg2./power, 'o-')
plot(wavenumbers(:,1), bg2, 'o-')
title('Background'); xlabel('Wavenumbers cm-1');

figure(1001)
[wavenumbers_sorted, I] = sort(wavenumbers(:,1));
pks2_sorted = pks2(I);
power = wavenumbers(:,2);
power_sorted = power(I);
pks2_sorted_power_corrected = pks2_sorted./power_sorted*mean(power_sorted);

plot(wavenumbers_sorted, pks2_sorted, '-o')
title('Peak height (power uncorrected)')
xlabel('Wavenumbers cm-1');

figure(1002)
plot(wavenumbers_sorted, pks2_sorted_power_corrected,'-o')
title('Peak height (power corrected)')
xlabel('Wavenumbers cm-1');





