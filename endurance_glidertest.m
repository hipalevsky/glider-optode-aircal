%% Analysis of glider air calibration test data from the Endurance array
addpath('G:\Shared drives\NSF_Irminger\Data_Files\Endurance_Glider_Test') %Processed glider data from Stuart Pierce, OSU
filename = '384-13_prototype_optode_data_v3.nc';

%% Load glider data
ce384.inair_oxygen_index = double(ncread(filename,'in-air_oxygen_index'));
info = ncinfo(filename);
for i=2:length(info.Variables)
    vn = info.Variables(i).Name;
    ce384.(vn) = ncread(filename, vn); 
end

%% Apply salinity and pressure corrections (oxygen provided is L1)
%Correcting assuming that internal salinity setting is 0
%Note: no lag correction applied in this processing
ce384.O2_corr = aaoptode_salpresscorr(ce384.oxygen, ce384.sea_water_temperature, ce384.salinity, ce384.pressure, 0);
ce384.O2sat_corr = ce384.saturation.*(1+ce384.pressure.*0.032./1000); %applies pressure but not salinity correction for saturation

%% Convert time
ce384.daten = datenum(1970,1,1,0,0,ce384.time);

%% Plot full depth data over entire deployments
plotstart = datenum(2023,6,15,0,0,0);
plotend = datenum(2023,8,11,0,0,0);

figure(1); clf
plot(ce384.daten, ce384.pressure, 'k.'); hold on;
scatter(ce384.daten, ce384.pressure, [], ce384.O2sat_corr,'filled'); colorbar; %caxis([48 122])
set(gca,'YDir','reverse'); 
xlim([plotstart plotend])
ylim([0 1000])
datetick('x',2,'keeplimits')
ylabel('Depth (m)')
title('Endurance Glider 384, oxygen saturation')

%% Calculate statistics for each air interval and keep only those with low stdev
    spacer = 10; %number of air measurements at beginning of each surfacing to cut
    tol = 1.5; %only keep air data if standard deviation of oxygen saturation measurements after cutting spacer is less than this value

%Number each of the surface air intervals
d = find(ce384.inair_oxygen_index == 1);
ce384.inair_profile = NaN*ce384.inair_oxygen_index;
ce384.inair_profile(d(1)) = 1;
for i = 2:length(d)
    if d(i) == d(i-1)+1
        ce384.inair_profile(d(i)) = ce384.inair_profile(d(i-1));
    else
        ce384.inair_profile(d(i)) = ce384.inair_profile(d(i-1)) + 1;
    end
end

intervals = unique(ce384.inair_profile(d));
ce384.airinterval_stats = NaN*ones(length(intervals), 5);
for i = 1:length(intervals)
    ind = find(ce384.inair_profile == intervals(i));
    ce384.airinterval_stats(i,1) = nanmean(ce384.saturation(ind(spacer + 1:end - spacer)));
    ce384.airinterval_stats(i,2) = nanstd(ce384.saturation(ind(spacer + 1:end - spacer)));
    ce384.airinterval_stats(i,3) = sum(~isnan(ce384.saturation(ind(spacer + 1:end - spacer))));
    ce384.airinterval_stats(i,4) = nanmean(ce384.daten(ind(spacer + 1:end - spacer)));
    if ce384.airinterval_stats(i,2) < tol & ce384.airinterval_stats(i,3) > 2 % count data as good if 3 or more measurements and std of O2sat < tolerance
        ce384.airinterval_stats(i,5) = 1; %put a "good" flag on usable data
    end
end
%% Visualize air interval data
figure(2); clf
histogram(ce384.saturation(d),[99:0.2:110]); hold on;
title('Histogram of all raw Endurance Glider 384 oxygen saturation data during air intervals')

figure(3); clf
plot(ce384.daten(d), ce384.saturation(d),'.', 'color', 'b'); hold on;
    ind = find(ce384.airinterval_stats(:,5) == 1);
plot(ce384.airinterval_stats(ind,4), ce384.airinterval_stats(ind,1), 'ro','markerfacecolor','r'); hold on;
datetick('x',2,'keeplimits')
legend('Glider 384, all','Glider 384, "good" interval mean','location','northeast')
ylabel('Raw oxygen saturation (%)')
title('Oxygen measurements during glider surface intervals measuring air, Endurance Array test of new mount')

%% Extract air and water measurements for each air interval

%%%% Constants used in calculations below
mbar2atm = 1013.25;
sec2day = 60*60*24;
    %Percentile of air measurement distribution to use
qntl = 0.32; %(value from Nicholson and Feen, 2017, p. 499)
    %Depths used to define near-surface oxygen measurements
surf_mindepth = 0.5;
surf_maxdepth = 10;
O2satmin = 20;
    %Time window for surface data: twin_beg sec < obs < twin_end sec
twin_beg = 90;
twin_end = 800;
    %Time window for corresponding water measurements
twin_wat = 15*60;

%Create a table to hold output
np = length(ce384.airinterval_stats); %number of air intervals
vars = {'air_daten','air_meas','air_corr','met_o2sat','ml_daten','ml_o2sat','ml_tem','ml_sal','nsurf'};
T = array2table(nan(np,length(vars)));
T.Properties.VariableNames = vars;

%Loop over all air intervals
for ii = 1:np   
    ind_air = find(ce384.inair_profile == ii);
    % find near-surface water measurements - use both dive and climb data if available (some just have dive)
    ind_wat = find(ce384.daten > (min(ce384.daten(ind_air)) - twin_wat/sec2day) & ce384.daten < (max(ce384.daten(ind_air)) + twin_wat/sec2day) & ...
        ce384.pressure < surf_maxdepth & ce384.pressure > surf_mindepth & ce384.saturation > O2satmin);
    
    % extract air oxygen measurements from designated time window during surface interval
    if length(ind_air) > spacer
        % time window for surface data: twin_beg sec < obs < twin_end sec
        t0 = ce384.daten(ind_air(1));
        ind_air_cut = intersect(ind_air, find((ce384.daten - t0) > twin_beg/sec2day & (ce384.daten - t0) < twin_end/sec2day));
        O2air = ce384.saturation(ind_air_cut);
        T.nsurf(ii) = length(ind_air_cut);
    else
        O2air = nan;
    end

    %Save air and water output value for each air interval
    T.ml_tem(ii) = nanmean(ce384.sea_water_temperature(ind_wat));
    T.ml_sal(ii) = nanmean(ce384.salinity(ind_wat));
    T.ml_o2sat(ii) = nanmedian(ce384.O2sat_corr(ind_wat));
    T.ml_daten(ii) = nanmean(ce384.daten(ind_wat));
    T.air_meas(ii) = quantile(O2air,qntl);
    T.air_daten(ii) = nanmean(ce384.daten(ind_air_cut));

    %Calculate O2 measurement expected in air based on meteorological data
    SVP_S = vpress(T.ml_sal(ii),T.ml_tem(ii)); % saturated water vapor pressure
    ph2o = nanmean(ce384.relative_humidity(ind_air_cut)).*SVP_S./100;
    T.met_o2sat(ii) = 100.*(nanmean(ce384.barometric_pressure/mbar2atm) - ph2o)./(1 - SVP_S);
end

%% Remove lines with missing data
d = ~isnan(T.air_meas+T.ml_o2sat) & T.air_meas > 0;
T(~d,:) = [];

%% Correct air measurements for surface water splashing
p = polyfit(T.ml_o2sat,T.air_meas,1);
T.air_corr = (T.air_meas-p(1).*T.ml_o2sat)./(1-p(1));

%% Calculate gain corrections

%This is a version not accounting for surface water splashing
px = polyfit(T.ml_daten-tref,T.met_o2sat./T.air_meas,2);
T.air_meas_corr = T.air_meas.*polyval(px,T.ml_daten-tref);
T.ml_o2sat_corr = T.ml_o2sat.*polyval(px,T.ml_daten-tref);

%Account for surface water splashing
px2 = polyfit(T.ml_daten-tref,T.met_o2sat./T.air_corr,1);
med_gain = median(T.met_o2sat(~isnan(T.met_o2sat))./T.air_corr(~isnan(T.met_o2sat)));

%% Figures
    %set figure options
ftsz = 14;
lnw = 1.5;
mrkr = 10;
gliderstring = 'Endurance G384';

%Plot time series of gain data with and without splash corrections
figure(4); clf
    ax1 = gca;
    hold all;
    ax1.FontSize = ftsz;
    dateplot = [floor(min(T.ml_daten)):ceil(max(T.ml_daten))];
plot(T.ml_daten, T.met_o2sat./T.air_meas, 'k.'); %not accounting for surface water splashing
    plot(dateplot, polyval(px,dateplot-tref),'k-','linewidth',lnw);
plot(T.ml_daten, T.met_o2sat./T.air_corr, 'b.'); %account for surface water splashing
    plot(dateplot, polyval(px2,dateplot-tref),'b-','linewidth',lnw);
    plot(dateplot, med_gain*ones(size(dateplot)),'b--','linewidth',lnw);
datetick('x',2,'keeplimits')
ylabel('Gain correction')
legend('Gain data w/ no splash correction','Gain quadratic fit w/ no splash correction',...
    'Gain data w/ splash correction','Gain linear fit w/ splash correction','Median gain w/ splash correction')
title([gliderstring ' gain corrections'])

%% Plot time series of glider data for gain calculations and corresponding MET data
figure(5); clf
    ax2 = gca;
    hold all;
    box on;
    cols = ax2.ColorOrder;
    ax2.FontSize = ftsz;
plot(T.ml_daten,T.ml_o2sat,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(1,:));
plot(T.air_daten,T.air_meas,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(3,:));
%plot(T.air_daten,T.air_corr,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(4,:));
plot(T.air_daten,med_gain.*T.air_corr,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color',cols(5,:));
plot(T.air_daten,T.met_o2sat,'.-','MarkerSize',mrkr,'LineWidth',lnw,'Color','k');
xlim([plotstart plotend])
title([gliderstring ' air calibration: gain = ' num2str(med_gain)]);
ylabel('Oxygen saturation (%)')
%ylim([98 106])
datetick('x',2,'keeplimits');
legend('\DeltaO_{2,w}^{meas}','\DeltaO_{2,a}^{meas}','\DeltaO_{2,a}^{splash & gain corr}','\DeltaO_{2}^{met}');

%% Plot relationship between glider surface water and air measurements
%(evidence of splashing)
figure(6); clf
    ax3 = gca;
    ax3.FontSize = ftsz;
    hold all;
    box on;
plot(T.ml_o2sat,T.air_meas,'.','MarkerSize',mrkr);
plot(ax3.XLim,p(1).*ax3.XLim+p(2),'-k','LineWidth',lnw);
xlabel('\DeltaO_{2,w}^{meas}');
ylabel('\DeltaO_{2,a}^{meas}');
    LM = fitlm(array2table([T.ml_o2sat,T.air_meas]),'linear');
title([gliderstring ' air vs. surface water, R^2 = ' num2str(LM.Rsquared.Adjusted,3) ', slope = ' num2str(table2array(LM.Coefficients(2,1)))]);

%% Relationship between MET and corrected glider air data
figure(7); clf
    ax4 = gca;
    ax4.FontSize = ftsz;
    hold all;
    box on;
plot(T.met_o2sat,T.air_corr,'.','MarkerSize',mrkr); hold on;
    ind = find(isnan(T.met_o2sat + T.air_corr) == 0);
    p2 = polyfit(T.met_o2sat(ind),T.air_corr(ind),1);
plot(ax4.XLim,p2(1).*ax4.XLim+p2(2),'-k','LineWidth',lnw);
xlabel('\DeltaO_{2}^{met}');
ylabel('\DeltaO_{2,a}^{splash & gain corr}');
    LM = fitlm(array2table([T.met_o2sat,T.air_corr]),'linear');
title([gliderstring ' corrected air vs. MET data, R^2 = ' num2str(LM.Rsquared.Adjusted,3)]);