% Figure 2 from Klavans et al. 2024
%
%
%
%


% Load files
time = ncread('Klavans_etal_PDO.nc','time');
PDO_OBS_NOAA = ncread('Klavans_etal_PDO.nc','PDO_OBS_NOAA');
PDO_OBS_NOAA_lp = ncread('Klavans_etal_PDO.nc','PDO_OBS_NOAA_lp');
PDO_EM = ncread('Klavans_etal_PDO.nc','PDO_EM');
PDO_EM_lp = ncread('Klavans_etal_PDO.nc','PDO_EM_lp');
PDO_E = ncread('Klavans_etal_PDO.nc','PDO_E');

%Internal only PDO
PDO_E = PDO_E(:,1:572);
for mm=1:length(PDO_E(1,:));
    PDO_int(:,mm) = PDO_E(:,mm) - PDO_EM;
end

%Forced ensembles
ens_size = 100;
cap = length(PDO_E(1,:));
for i=1:572;
    mem = ceil(cap*rand(ens_size,1));
    PDO_EM_rand(:,i) = mean(PDO_E(:,mem),2);
end

edges = -1:0.05:1;
for i=1:length(PDO_E(1,:));
    PDO_tmp = jk_filter(PDO_E(:,i),10);
    rhist(i) = corr(PDO_tmp,PDO_OBS_NOAA_lp);
    PDO_tmp_uf = jk_filter(PDO_int(:,i),10);
    rhist_uf(i) = corr(PDO_tmp_uf,PDO_OBS_NOAA_lp);
    PDO_tmp_f = jk_filter(PDO_EM_rand(:,i),10);
    rhist_f(i) = corr(PDO_tmp_f,PDO_OBS_NOAA_lp);
end
rhist_EM = corr(jk_filter(PDO_EM,10),PDO_OBS_NOAA_lp);

rng(823,'twister');
btstrp = 1000; %change to 10,000 for published results
indx = ceil(572*rand(1,1));

for i=1:500;
    for j=1:btstrp;
        mem = ceil(572*rand(i,1));
        EM_rand = mean(PDO_E(:,mem),2);
        rhist_obs(i,j) = corr(jk_filter(EM_rand,10),PDO_OBS_NOAA_lp);
    end
end
rhist_obs_sort = sort(rhist_obs,2);



%% Plot
fig1 = figure(1);
set(fig1,'color','w','units','inches');
set(fig1,'pos',[3 3 11 8]);

subplot(2,3,4);
hold on
p2 = plot([rhist_EM rhist_EM],[0 50],'color',[0.2 0.2 0.2],'LineWidth',5);
p1 = histogram(rhist_uf,edges,'FaceColor',[.35 .6 .9],'Edgecolor',[.35 .6 .9],'Normalization','probability');
xlim([-1 1]);
ylim([0 .31]);
grid on; box on;
xlabel('Correlation coefficient');
ylabel('Probability');
set(gca,'YTick',0:0.1:0.3);
text(-.98,.325,'d.','FontWeight','bold','FontSize',13);
legend([p2,p1],{'Forcing','Internal'},'Location','NorthWest');
legend boxoff
set(gca,'FontSize',13);

subplot(2,3,5);
hold on
p4 = plot([rhist_EM rhist_EM],[0 50],'color',[0.2 0.2 0.2],'LineWidth',5);
p3 = histogram(rhist,edges,'FaceColor',[.35 .6 .9],'Edgecolor',[.35 .6 .9],'Normalization','probability');
xlim([-1 1]);
ylim([0 .31]);
grid on; box on;
xlabel('Correlation coefficient');
%ylabel('Probability');
set(gca,'YTick',0:0.1:0.3);
set(gca,'YTickLabel',{'' '' '' ''});
text(-.98,.325,'e.','FontWeight','bold','FontSize',13);
legend([p4,p3],{'Forcing','Int. + Forcing'},'Location','NorthWest');
legend boxoff
set(gca,'FontSize',13);

subplot(2,3,6);
hold on
p4 = plot([rhist_EM rhist_EM],[0 200],'color',[0.2 0.2 0.2],'LineWidth',5);
p3 = histogram(rhist_f,edges,'FaceColor',[.35 .6 .9],'Edgecolor',[.35 .6 .9],'Normalization','probability');
xlim([-1 1]);
ylim([0 .31]);
grid on; box on;
xlabel('Correlation coefficient');
%ylabel('Probability');
set(gca,'YTick',0:0.1:0.3);
set(gca,'YTickLabel',{'' '' '' ''});
text(-.98,.325,'f.','FontWeight','bold','FontSize',13);
legend([p4,p3],{'Forcing','Forcing'},'Location','NorthWest');
legend boxoff
set(gca,'FontSize',13);

% SN plot

subplot(2,3,1);
hold on
scatter(1:500,100*mean(rhist_obs,2).^2,'MarkerFaceColor',[.35 .6 .9]);
ciplot(100*rhist_obs_sort(:,round(0.975*btstrp)).^2,100*rhist_obs_sort(:,round(0.025*btstrp)).^2,1:500,[.35 .6 .9]);
ylabel('Explained variance (%)','FontSize',13);
xlabel('Ensemble size','FontSize',13)
ylim([-5 70]);
xlim([1 500]);
set(gca,'Ytick',0:25:100);
text(5,73,'a.','FontWeight','bold','FontSize',13);
grid on; box on;
set(gca,'FontSize',13);


%% 1850 ensemble
clearvars -except fig1 rhist_EM rhist_obs

time = ncread('Klavans_etal_PDO_1870ens.nc','time');
PDO_OBS_NOAA = ncread('Klavans_etal_PDO_1870ens.nc','PDO_OBS_NOAA');
PDO_OBS_NOAA_lp = ncread('Klavans_etal_PDO_1870ens.nc','PDO_OBS_NOAA_lp');
PDO_GrandEM = ncread('Klavans_etal_PDO_1870ens.nc','PDO_GrandEM');
PDO_GrandEM_lp = ncread('Klavans_etal_PDO_1870ens.nc','PDO_GrandEM_lp');
PDO_GrandE = ncread('Klavans_etal_PDO_1870ens.nc','PDO_GrandE');
color_list = [.35 .7 .9;...
              0 .45 .7;...
              .35 .6 .9;...
              .4 .8 0];

%%1920 test
for kk=1:1000;
    tmp = jk_filter(mean(PDO_GrandE(71:end,ceil(381*rand(40,1))),2),10);
    r1920(kk) = corr(tmp,PDO_OBS_NOAA_lp(51:end));
end

%Filtered correlations
r_1870_lp = corr(PDO_GrandEM_lp(21:end),PDO_OBS_NOAA_lp(1:end));
r_1870_1950_lp = corr(PDO_GrandEM_lp(21:101),PDO_OBS_NOAA_lp(1:81));
r_1900_lp = corr(PDO_GrandEM_lp(51:end),PDO_OBS_NOAA_lp(31:end));
r_1925_lp = corr(PDO_GrandEM_lp(76:end),PDO_OBS_NOAA_lp(56:end));
 
bar_val_lp = 100*[r_1870_lp r_1870_1950_lp rhist_EM].^2;

rng(823,'twister');
btstrp = 1000; 

%1870-2014
for j=1:btstrp;
    mem = ceil(381*rand(100,1));
    EM_rand = mean(PDO_GrandE(21:end,mem),2);
    rhist_full(j) = corr(jk_filter(EM_rand,10),PDO_OBS_NOAA_lp(1:end));
end
rhist_full_sort = sort(rhist_full,2);

%1870-1950
for j=1:btstrp;
    mem = ceil(381*rand(100,1));
    EM_rand = mean(PDO_GrandE(21:101,mem),2);
    rhist_early(j) = corr(jk_filter(EM_rand,10),PDO_OBS_NOAA_lp(1:81));
end
rhist_early_sort = sort(rhist_early,2);

%1870-1950
for j=1:btstrp;
    mem = ceil(381*rand(100,1));
    EM_rand = mean(PDO_GrandE(102:end,mem),2);
    rhist_late(j) = corr(jk_filter(EM_rand,10),PDO_OBS_NOAA_lp(82:end));
end
rhist_late_sort = sort(rhist_late,2);

err_val_low = [rhist_full_sort(25)-r_1870_lp rhist_early_sort(25)-r_1870_1950_lp rhist_late_sort(25)-rhist_EM];
err_val_high = [rhist_full_sort(975)-r_1870_lp rhist_early_sort(975)-r_1870_1950_lp rhist_late_sort(975)-rhist_EM];
if r_1870_lp+err_val_low(1) < 0
    err_val_low(1) = r_1870_lp;
end
if r_1870_1950_lp+err_val_low(2) < 0
    err_val_low(2) = r_1870_1950_lp;
end
if rhist_EM+err_val_low(3) < 0
    err_val_low(3) = rhist_EM;
end

err_val_low = 100*err_val_low.^2;
err_val_high = 100*err_val_high.^2;
    
subplot(2,3,2);
hold on
b1 = bar(1:3,bar_val_lp,'FaceColor','flat');
for k = 1:size(bar_val_lp,2)
    b1.CData(k,:) = color_list(k,:);
end

er = errorbar(1:3,bar_val_lp,err_val_low,err_val_high);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

tickLabels = {'1870-2014','1870-1950','1950 - 2014'};

xlim([.5 3.5]);
ylim([-5 70]);
ax = gca(); 
ax.XTick = 1:3; 
ax.XTickLabel = tickLabels; 
set(gca,'Ytick',0:25:100);
set(gca,'FontSize',13);
%ylabel('Explained variance (%)','FontSize',13);
text(0.55,73,'b.','FontWeight','bold','FontSize',13);
grid on
box on
xtickangle(45)
%pos2 = get(s2,'pos');
%set(s2,'pos',[pos2(1) pos2(2)+0.01 pos2(3)-0.102 pos2(4)]);

%% Physics ensembles
ens_list = {'CMIP5','CMIP6','EMISSIONS','CONCENTRATIONS','INTERACTIVE','noINTERACTIVE'};
title_list = {'CMIP5','CMIP6','Emissions','Concentrations','Interactive','Not interactive'};
title_list2 = {'270','302','460','112','442','130'};

time = ncread('Klavans_etal_PDO.nc','time');
PDO_OBS_NOAA = ncread('Klavans_etal_PDO.nc','PDO_OBS_NOAA');
PDO_OBS_NOAA_lp = ncread('Klavans_etal_PDO.nc','PDO_OBS_NOAA_lp');


for ee=1:length(ens_list)

    ens_name = ens_list{ee};

    PDO_EM = ncread(['Klavans_etal_PDO_',ens_name,'.nc'],'PDO_EM');
    PDO_EM_lp = ncread(['Klavans_etal_PDO_',ens_name,'.nc'],'PDO_EM_lp');
    PDO_E = ncread(['Klavans_etal_PDO_',ens_name,'.nc'],'PDO_E');
    
    bar_val_lp(ee) = 100*corr(PDO_EM_lp,PDO_OBS_NOAA_lp).^2;

    for j=1:btstrp;
        mem = ceil(length(PDO_E(1,:))*rand(100,1));
        EM_rand = mean(PDO_E(1:end,mem),2);
        rhist_ens(ee,j) = corr(jk_filter(EM_rand,10),PDO_OBS_NOAA_lp(1:end));
    end
    rhist_ens_sort = sort(rhist_ens,2);

    eval_l_tmp = rhist_ens_sort(25)-corr(PDO_EM_lp,PDO_OBS_NOAA_lp);
    eval_h_tmp = rhist_ens_sort(975)-corr(PDO_EM_lp,PDO_OBS_NOAA_lp);
    eval_l_tmp(corr(PDO_EM_lp,PDO_OBS_NOAA_lp)+eval_l_tmp<0) = corr(PDO_EM_lp,PDO_OBS_NOAA_lp);

    er_val_low_ens(ee) = 100*eval_l_tmp^2;
    er_val_high_ens(ee) = 100*eval_h_tmp^2;

end

subplot(2,3,3);
hold on
b2 = bar(1:6,bar_val_lp,'FaceColor','flat');
for k = 1:size(bar_val_lp,2)
    b2.CData(k,:) = [.35 .6 .9];
end
hT=[];              % placeholder for text object handles
for i=1:6  % iterate over number of bar objects
  hT=[hT,text(b2.XData(i),double(b2.YData(i))+1,title_list2(:,i), ...
          'VerticalAlignment','bottom','horizontalalign','center')];
end

labelArray = [title_list];
tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));

er = errorbar(1:6,bar_val_lp,er_val_low_ens,er_val_high_ens);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';

xlim([.5 6.5]);
ylim([-5 70]);
ax = gca(); 
ax.XTick = 1:6; 
ax.XTickLabel = labelArray; 
set(gca,'Ytick',0:25:100);
set(gca,'FontSize',13);
%ylabel('Explained variance (%)','FontSize',13);
text(0.55,73,'c.','FontWeight','bold','FontSize',13);
grid on
box on
xtickangle(45)



