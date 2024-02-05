% Figure 3 from Klavans et al. PDO paper
%
%
%
%


% set up plot
fig = figure(1);
set(fig, 'color', 'w','units','inches'); 
set(fig,'pos',[3 3 11 8.5]);

color_list_bar = [.95 .9 .25
              .9 .6 0;...
              .8 .4 0;...
              .8 .6 .7;...
              .4 .8 0;...
              .4 .5 .8];
tickLabels = {'All-forcings','Aerosol-only','GHG-only','Natural-only'};

% Load Files
time = ncread('Klavans_etal_PDO.nc','time');
PDO_OBS_NOAA = ncread('Klavans_etal_PDO.nc','PDO_OBS_NOAA');
PDO_OBS_NOAA_lp = ncread('Klavans_etal_PDO.nc','PDO_OBS_NOAA_lp');
PDO_EM = ncread('Klavans_etal_PDO.nc','PDO_EM');
PDO_EM_lp = ncread('Klavans_etal_PDO.nc','PDO_EM_lp');

% all forcing
r_full_lp = corr(PDO_EM_lp,PDO_OBS_NOAA_lp);
r_pre80_lp = corr(PDO_EM_lp(1:40),PDO_OBS_NOAA_lp(1:40));
r_post80_lp = corr(PDO_EM_lp(41:end),PDO_OBS_NOAA_lp(41:end));

bar_val_lp = 100*[r_pre80_lp r_post80_lp r_full_lp].^2;
    
s5 = subplot(2,1,2);
hold on
b1 = bar(1,bar_val_lp,'FaceColor','flat');
for k = 1:size(bar_val_lp,2)
    b1(k).CData = color_list_bar(k,:);
end

xlim([.5 4.5]);
ylim([-5 70]);
ax = gca(); 
ax.XTick = 1:4; 
ax.XTickLabel = tickLabels; 
set(gca,'Ytick',0:25:100);
set(gca,'FontSize',13);
ylabel('Explained variance (%)','FontSize',13);
text(0.55,65,'d.','FontWeight','bold','FontSize',13);
grid on
box on
xtickangle(45)
pos5 = get(s5,'pos');
%set(s2,'pos',[pos2(1) pos2(2)+0.01 pos2(3)-0.136 pos2(4)]);

clearvars -except fig color_list_bar PDO_OBS_NOAA PDO_OBS_NOAA_lp

% Single forcing timeseries
scenario_list = {'hist-aer','hist-GHG','hist-nat'};
label_list = {'Aerosol-only','GHG-only','Natural-only'};
letter_list = {'a.','b.','c.'};
color_list = [0 .45 .7;...
              .6 .9 0;...
              .8 .6 .7;...
              .4 .8 0;...
              0 .6 .5;...
              .95 .9 .25;...
              .35 .7 .9];
PDO_OBS_NOAA = PDO_OBS_NOAA(1:end-1);
PDO_OBS_NOAA_lp = PDO_OBS_NOAA_lp(1:end-1);

hold on
for iii=1:length(scenario_list);
        scenario = scenario_list{iii};

        PDO_EM = ncread(['Klavans_etal_PDO_',scenario,'.nc'],'PDO_EM');
        PDO_EM_lp = ncread(['Klavans_etal_PDO_',scenario,'.nc'],'PDO_EM_lp');
        PDO_E = ncread(['Klavans_etal_PDO_',scenario,'.nc'],'PDO_E');
        for i=1:75;
            PDO_E_lp(:,i) = jk_filter(PDO_E(:,i),10);
        end

        load(['pdo_singForcing_',scenario,'_1950_2013_NEW.mat']);
        PDO_EM_norm = (PDO_EM-mean(PDO_EM))/std(PDO_EM);
        PDO_E_sort = sort(PDO_E,2);
        
        subplot(2,3,iii);
        plot([1950 2013],[0 0],'k');
        hold on
        plot(1950:2013,PDO_OBS_NOAA,'k','linewidth',2)
        plot(1950:2013,PDO_OBS_NOAA_lp,'k','linewidth',5)
        plot(1950:2013,PDO_EM,'color',[.5 .5 .5],'linewidth',2);
        plot(1950:2013,PDO_EM_norm,'color',color_list(iii,:),'linewidth',2); 
        plot(1950:2013,jk_filter(PDO_EM_norm,10),'color',color_list(iii,:),'linewidth',5);
        box on;grid on
        if iii == 1
            ylabel('PDO Index');
        end
        xlabel('Year');
        set(gca,'FontSize',13);
        ylim([-4 4]);
        xlim([1950 2013]);
        text(1951,3.6,[letter_list{iii},' ',label_list{iii}],'FontWeight','bold','FontSize',14);
        
        r_full_lp = corr(jk_filter(PDO_EM,10),jk_filter(PDO_OBS_NOAA,10));
        [mu,ciLo_full,ciHi_full] = ens_bootstrap(PDO_OBS_NOAA_lp,PDO_E_lp,10000);
        r_pre80_lp = corr(jk_filter(PDO_EM(1:40),10),jk_filter(PDO_OBS_NOAA(1:40),10));
        [mu,ciLo_pre80,ciHi_pre80] = ens_bootstrap(PDO_OBS_NOAA_lp(1:40),PDO_E_lp(1:40,:),10000);
        r_post80_lp = corr(jk_filter(PDO_EM(41:end),10),jk_filter(PDO_OBS_NOAA(41:end),10));
        [mu,ciLo_post80,ciHi_post80] = ens_bootstrap(PDO_OBS_NOAA_lp(41:end),PDO_E_lp(41:end,:),10000);
   
        bar_val_lp = 100*[r_pre80_lp r_post80_lp r_full_lp].^2
        %ci_val_lo = 100*[ciLo_full ciLo_pre80 ciLo_post80].^2
        ci_val_hi = 100*[ciHi_full ciHi_pre80 ciHi_post80].^2
        %bar_val_lp(bar_val_lp<0) = NaN;
    
        s5 = subplot(2,1,2);
        hold on
        b1 = bar(iii+1,bar_val_lp,'FaceColor','flat');
        for k = 1:size(bar_val_lp,2)
            b1(k).CData = color_list_bar(k,:);
        end

        clear('PDO_E','PDO_EM','r_PDO','PDO_INDEX','bar_val_lp','PDO_EM_lp','PDO_INDEX_OBS_lp');
        
end

legend(b1,{'PDO (1950 - 1989)','PDO (1990 - 2014)',...
    'PDO (1950 - 2014)'},...
    'Location','NorthEast');
legend boxoff

