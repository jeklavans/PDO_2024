% Figure 1 from Klavans et al. (2024)
%
%
%
%

% Load Files
time = ncread('Klavans_etal_PDO.nc','time');
PDO_OBS_NOAA = ncread('Klavans_etal_PDO.nc','PDO_OBS_NOAA');
PDO_OBS_NOAA_lp = ncread('Klavans_etal_PDO.nc','PDO_OBS_NOAA_lp');
PDO_EM = ncread('Klavans_etal_PDO.nc','PDO_EM');
PDO_EM_lp = ncread('Klavans_etal_PDO.nc','PDO_EM_lp');

latSST_OBS = ncread('Klavans_etal_PDO.nc','latSST_OBS');
lonSST_OBS = ncread('Klavans_etal_PDO.nc','lonSST_OBS');
latPSL_OBS = ncread('Klavans_etal_PDO.nc','latPSL_OBS');
lonPSL_OBS = ncread('Klavans_etal_PDO.nc','lonPSL_OBS');
beta_SST_PDO_OBS = ncread('Klavans_etal_PDO.nc','beta_SST_PDO_OBS');
beta_PSL_PDO_OBS = ncread('Klavans_etal_PDO.nc','beta_PSL_PDO_OBS');

latSST = ncread('Klavans_etal_PDO.nc','latSST');
lonSST = ncread('Klavans_etal_PDO.nc','lonSST');
latPSL = ncread('Klavans_etal_PDO.nc','latPSL');
lonPSL = ncread('Klavans_etal_PDO.nc','lonPSL');
beta_SST_PDO = ncread('Klavans_etal_PDO.nc','beta_SST_PDO');
beta_PSL_PDO = ncread('Klavans_etal_PDO.nc','beta_PSL_PDO');

% Set up plot
fig = figure(1);
set(fig, 'color', 'w','units','inches'); 
set(fig,'pos',[3 3 10 8]);

box1 = [75 -10 100 280];
rb_short = [103,0,31;...
    178,24,43;...
    214,96,77;...
    244,165,130;...
    253,219,199;...
    247,247,247;...
    209,229,240;...
    146,197,222;...
    67,147,195;...
    33,102,172;...
    5,48,97];
rb_short = rb_short./256;
rb_short = flipud(rb_short);

% Plot timeseries
subplot(2,1,1);
hold on
plot([1950 2014],[0 0],'k');
plot(time,PDO_OBS_NOAA,'k','LineWidth',2);
p1 = plot(time,PDO_OBS_NOAA_lp,'k','LineWidth',5);
ylim([-2 2]);
xlim([1950 2014]);
grid on; box on;
text(1950.5,1.85,'a.','FontSize',13,'FontWeight','bold')
ylabel('PDO index','FontSize',13);

%normalize for viz
PDO_EM_norm = (PDO_EM-mean(PDO_EM))/std(PDO_EM);
PDO_EM_lp_norm = (PDO_EM_lp-mean(PDO_EM_lp))/std(PDO_EM_lp);

s1 = subplot(2,1,1);
p2a=plot(time,PDO_EM,'color',[5,113,176]./256,'LineWidth',3);
plot(time,PDO_EM_norm,'color',[.35 .7 .9 0.5],'LineWidth',2);
p2 = plot(time,PDO_EM_lp_norm,'color',[.35 .7 .9],'LineWidth',5);
[r,p] = corrcoef(jk_filter(PDO_EM,10),jk_filter(PDO_OBS_NOAA,10))
text(1950.5,-1.82,['R^2 = ',num2str(100*r(1,2)^2,2),'%'],'FontSize',13);
pos1 = get(s1,'pos');
xlabel('Year');
set(gca,'FontSize',13);

lgd = legend([p1,p2],{'Obs.','Models'},'Location','NorthEast');
set(lgd,'FontSize',13);
legend boxoff

% PDO Patterns

%% SST Patterns

%observed PDO pattern
latSST_OBS = double(latSST_OBS);
lonSST_OBS = double(lonSST_OBS);
[mlat,mlon] = meshgrid(latSST_OBS,lonSST_OBS); 
mlat = [mlat; mlat(1,:)];
mlon = [mlon; mlon(1,:)];
latPSL_OBS = double(latPSL_OBS); 
lonPSL_OBS = double(lonPSL_OBS);

s3 = subplot(2,2,3);
SetMap_NAtlantic;
colormap(gca,rb_short);
pcolorm(mlat,mlon,-1*beta_SST_PDO_OBS,'LineStyle','None');
hold on
[C,h] = contourm(latPSL_OBS,lonPSL_OBS,-1*beta_PSL_PDO_OBS',0,'k','LineWidth',2); 
[C,h] = contourm(latPSL_OBS,lonPSL_OBS,-1*beta_PSL_PDO_OBS',-3:.5:-.1,'color',[153,142,195]./256,'LineWidth',2);
caxis([-0.5 0.5]);
SetMap_NAtlantic_Gray;
ax = axesm('eqdcylin','Grid','off',...
    'Frame','on',...
    'ParallelLabel','on',...
    'MeridianLabel','on',...
    'MLabelParallel','south',...
    'MLabelLocation',[-270 -225 -180 -135 -90 -45 0],...
    'PLabelLocation',[-60 -30 0 30 60],...
    'MapLatLimit',[box1(1) box1(2)],...
    'MapLonLimit',[box1(3) box1(4)]);
linem([36 36 31 31 36],[140 165 165 140 140],'k','LineWidth',2);

h = mlabel('on');
for i=1:length(h)
  labelStrings(i)=string(h(i).String{2});
  value = str2double(extractBetween(labelStrings(i),2,"^"));
  value = wrapTo360(value);
  if value > 180
    value = 360 - value;    
    labelStrings(i) = strcat(" ",num2str(value),"^{\circ} E");
    h(i).String{2} = labelStrings(i);
  end
end

y1 = title('b. Observations','FontSize',13,'FontName','Helvetica');
py1 = get(y1,'pos');
set(y1,'pos',[py1(1)-0.75,py1(2),py1(3)]);
set(gca,'FontName','Helvetica','FontSize',13);

s3p = get(s3,'pos');
set(s3,'pos',[s3p(1) s3p(2)+0.07 s3p(3) s3p(4)]);


%% forced patterns
latSST = double(latSST);
lonSST = double(lonSST);
[mlat,mlon] = meshgrid(latSST,lonSST); 
mlat = [mlat; mlat(1,:)];
mlon = [mlon; mlon(1,:)];

s4 = subplot(2,2,4);
SetMap_NAtlantic;
colormap(gca,rb_short);
pcolorm(mlat,mlon,beta_SST_PDO,'LineStyle','None');
hold on
latPSL = ncread('Klavans_etal_PDO.nc','latPSL');
[C,h] = contourm(latPSL,lonPSL,beta_PSL_PDO',0,'k','LineWidth',2);
[C,h] = contourm(latPSL,lonPSL,beta_PSL_PDO',-3:.5:-0.1,'color',[153,142,195]./256,'LineWidth',2);
caxis([-0.5 0.5]); 
SetMap_NAtlantic_Gray;
ax = axesm('eqdcylin','Grid','off',...
    'Frame','on',...
    'ParallelLabel','on',...
    'MeridianLabel','on',...
    'MLabelParallel','south',...
    'MLabelLocation',[-270 -225 -180 -135 -90 -45 0],...
    'PLabelLocation',[-60 -30 0 30 60],...
    'MapLatLimit',[box1(1) box1(2)],...
    'MapLonLimit',[box1(3) box1(4)]);

h = mlabel('on');
for i=1:length(h)
  labelStrings(i)=string(h(i).String{2});
  value = str2double(extractBetween(labelStrings(i),2,"^"));
  value = wrapTo360(value);
  if value > 180
    value = 360 - value;    
    labelStrings(i) = strcat(" ",num2str(value),"^{\circ} E");
    h(i).String{2} = labelStrings(i);
  end
end

y1 = title('c. Models','FontSize',13,'FontName','Helvetica');
py1 = get(y1,'pos');
set(y1,'pos',[py1(1)-1,py1(2),py1(3)]);
set(gca,'FontName','Helvetica');

s4p = get(s4,'pos');
set(s4,'pos',[s4p(1) s4p(2)+0.07 s4p(3) s4p(4)]);

cbar = colorbar('Location','SouthOutside');
pbar = get(cbar,'pos');
set(cbar,'pos',[pbar(1)-0.37 pbar(2)-0.08 pbar(3)+0.3 0.01]);
set(cbar,'Ticks',[-.5 0 .5],'Ticklength',0);
xlabel(cbar,'Regression coefficient (^{\circ}C^ unit PDO^{-1})','FontSize',13,'FontName','Helvetica');
set(gca,'FontName','Helvetica','FontSize',13);



