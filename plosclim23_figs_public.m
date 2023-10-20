function plosclim23_figs_public
% =========================================================================
% Plot figures from Jacox et al. (2023). See paper for interpretation.
%
% Jacox, M.G., M. Pozo Buil, S. Brodie, M.A. Alexander, D.J. Amaya, S.J.
% Bograd, C.A. Edwards, J. Fiechter, E.L. Hazen, G. Hervieux, D. Tommasi
% (2023), Downscaled seasonal forecasts for the California Current System:
% Skill assessment and prospects for living marine resource applications,
% PLOS Climate, doi:10.1371/journal.pclm.0000245.
%
% Each function loads in the appropriate data for that figure, e.g.
% fig1.mat. Download the data files and place them in the same folder as
% this script.
%
% M. Jacox
% October 2023
% =========================================================================

% Uncomment figures to plot
% ==================
% fig1 % Sample SST forecast output
% fig2  % Forecast skill summary for all variables, lead times, and metrics
% fig3 % Skill maps for all variables
% fig4 % SST skill maps
% fig5 % SSH skill maps
% fig6('a') % Bottom temp skill maps 
% %   Panels are plotted individually
% %   Input panel letter 'a' 'b' 'c' or 'd', e.g., fig6('a')
% %   a: Jan initialization, 1.5-month lead
% %   b: Jan initialization, 6.5-month lead
% %   c: Jul initialization, 1.5-month lead
% %   d: Jul initialization, 6.5-month lead
% fig7 % Evaluate SST forecasts against shore stations
% fig8 % Evaluate SST forecasts against shore stations
% fig9 % Compare forecast ensemble spread to observed variability
% fig10 % Compare skill of individual ensemble members for all variables
% fig11 % Forecast skill summary for all variables and metrics, int vs. bc
% fig12 % Skill for specific LMR applications, including reduced versions with only skillful variables
% fig13 % Compare skill above persistence for sst, ssh, and wind
% fig14 % SSH and SST skill maps with verification from sensitivity runs

% ==========================
function fig1
% ==========================
% Sample forecast output

% Load input file
f_in = 'fig1.mat';
load(f_in)

% Plotting info
c = [1 .6 0;0 .7 0]; % Colors for forecast plots
pos = [.05 .1 .3 .75;... % Position of panels in plot
       .45 .5 .48 .35;...
       .45 .1 .48 .35];
cbar_pos = [pos(1,1) pos(1,2)+pos(1,4)+.02 pos(1,3) .02];
dist_cont = [75 300]; % offshore distances dividing regions
lat_cont = [30.5 34.5 40.5 47.5]; % latitudes dividing regions
init_mon = 1; % January initializations

% ========================
% Plot SST
figure,set(gcf,'color','w','position',[100 100 900 450])
hold on

% PLOT MEAN SST MAP
subplot(221)
hold on
pcolor(lon-.05,lat-.05,sst_obs_mean)
shading flat
set(gca,'fontsize',10,'tickdir','out')
axis([-134 -115.5 30 48])
box on
set(gca,'color',[.8 .8 .8])
caxis([10 20])
cmocean('thermal',20)
cb = colorbar;
set(cb,'fontsize',12,'tickdir','out','orientation','horizontal')
xlabel(cb,['SST (' char(176) 'C)'])
set(gca,'position',pos(1,:))
set(cb,'position',cbar_pos)

% Plot coastline and contours
for ii = 1:length(dist_cont)
    for jj = 1:length(lat_cont)-1  
        mask = zeros(size(dist_coast));
        mask(dist_coast>0 & dist_coast<=dist_cont(ii) & lat>=lat_cont(jj) & lat<=lat_cont(jj+1)) = 1;
        mask(isnan(dist_coast)) = nan;
        contour(lon,lat,mask,[.99 .99],'k','linewidth',1)
    end
end
plot(coastline(:,1),coastline(:,2),'k.','markersize',4)

% Plot observed time series
subplot(222),hold on
h(3) = plot(decyear_obs,sst_obs_reg,'k','linewidth',3);
subplot(224),hold on
plot(decyear_obs,ssta_obs_reg,'k','linewidth',3)

% Plot forecasts
for ii = 1:length(year_forc)
    
    % Forecast time vector
    time_init = year_forc(ii)+init_mon/12-1/24;
    decyear_forc = time_init:1/12:time_init+11/12;
    
    % PLOT SST
    subplot(222)
    
    % Plot individual ensemble members
    fill([decyear_forc fliplr(decyear_forc)],[min(squeeze(sst_int(:,:,ii))) fliplr(max(squeeze(sst_int(:,:,ii))))],c(1,:),'edgecolor','none','facealpha',0.3)
    fill([decyear_forc fliplr(decyear_forc)],[min(squeeze(sst_bc(:,:,ii))) fliplr(max(squeeze(sst_bc(:,:,ii))))],c(2,:),'edgecolor','none','facealpha',0.3)
    
    % Plot ensemble means
    h(1) = plot(decyear_forc,squeeze(mean(sst_int(:,:,ii))),'linewidth',2,'color',c(1,:));
    h(2) = plot(decyear_forc,squeeze(mean(sst_bc(:,:,ii))),'linewidth',2,'color',c(2,:));
    
    % Plot forecast starts
    plot(decyear_forc(1),mean(sst_int(:,1,ii)),'o','color',c(1,:),'markerfacecolor',c(1,:),'markersize',7)
    plot(decyear_forc(1),mean(sst_bc(:,1,ii)),'o','color',c(2,:),'markerfacecolor',c(2,:),'markersize',7)
    
    % Format
    if ii==length(year_forc)
        axis([2000 2010 10 18])
        box on
        grid on
        set(gca,'tickdir','out','fontsize',12)
        ylabel(['SST (' char(176) 'C)'])
        set(gca,'xticklabel',[])
        set(gca,'position',pos(2,:))
        hl = legend(h,'CanCM4-ROMS','CanCM4-ROMS-BC','Observed');
        set(hl,'fontsize',12,'orientation','horizontal')
        
        % Center legend
        hl_width = hl.Position(3);
        xspace = pos(2,3) - hl_width; % Used to center legend
        set(hl,'fontsize',12,'orientation','horizontal','position',[pos(2,1)+xspace/2 pos(2,2)+pos(2,4)+.02 hl_width .05])
    end
    
    % PLOT SSTa
    subplot(224)
    
    % Plot individual ensemble members
    fill([decyear_forc fliplr(decyear_forc)],[min(squeeze(ssta_int(:,:,ii))) fliplr(max(squeeze(ssta_int(:,:,ii))))],c(1,:),'edgecolor','none','facealpha',0.3)
    fill([decyear_forc fliplr(decyear_forc)],[min(squeeze(ssta_bc(:,:,ii))) fliplr(max(squeeze(ssta_bc(:,:,ii))))],c(2,:),'edgecolor','none','facealpha',0.3)
    
    % Plot ensemble means
    plot(decyear_forc,squeeze(mean(ssta_int(:,:,ii))),'linewidth',2,'color',c(1,:))
    plot(decyear_forc,squeeze(mean(ssta_bc(:,:,ii))),'linewidth',2,'color',c(2,:))
    
    % Plot forecast starts
    plot(decyear_forc(1),mean(ssta_int(:,1,ii)),'o','color',c(1,:),'markerfacecolor',c(1,:),'markersize',7)
    plot(decyear_forc(1),mean(ssta_bc(:,1,ii)),'o','color',c(2,:),'markerfacecolor',c(2,:),'markersize',7)
    
    % Format
    if ii==length(year_forc)
        axis([2000 2010 -2 2])
        box on
        grid on
        set(gca,'tickdir','out','fontsize',12)
        ylabel(['SST anomaly (' char(176) 'C)'])
        set(gca,'position',pos(3,:))
    end
end

% ==========================
function fig2
% ==========================
% Skill summary for all variables and metrics

% Load input file
f_in = 'fig2.mat';
load(f_in)

% Variables names
vars = {'curl' 'svstr' 'sustr' 'sv' 'su' 'mld'...
        'bf' 'bt' 'ssh' 'sst'};
    
% Plotting info
xpos = [.08 .5 .08 .5 .08 .5 .08 .5];
ypos = [.75 .75 .53 .53 .31 .31 .09 .09];
width = .37;
height = .19;

% Open figure
figure,set(gcf,'color','w','position',[100 100 475 750])

% Plot
[ni,nv,nl] = size(acc);

for ii = 1:ni
    
    subplot(4,ni,ii),hold on
    varplot = nan(nv+1,nl+1); % Pad for plotting with pcolor
    sig_zero_plot = nan(nv+1,nl+1);
    sig_pers_plot = nan(nv+1,nl+1);
    varplot(1:end-1,1:end-1) = squeeze(acc(ii,:,:));
    sig_zero_plot(1:end-1,1:end-1) = squeeze(acc_sig_zero(ii,:,:));
    sig_pers_plot(1:end-1,1:end-1) = squeeze(acc_sig_pers(ii,:,:));
    varplot(varplot<0) = nan;
    pcolor(lead,varnum,varplot)
    set(gca,'color','w')
    caxis([0 1])
    cmocean('amp',10)
    ind = find(sig_zero_plot==1);
    plot(lead(ind)+.5,varnum(ind)+.5,'o','color',[.6 .6 .6],'markerfacecolor',[.6 .6 .6],'markersize',5)
    ind = find(sig_zero_plot==1 & sig_pers_plot==0);
    plot(lead(ind)+.5,varnum(ind)+.5,'wo','markerfacecolor','w','markersize',5)
    
    subplot(4,ni,ii+2),hold on
    varplot = nan(nv+1,nl+1);
    varplot(1:end-1,1:end-1) = squeeze(fa(ii,:,:));
    varplot(varplot<0.56) = nan;
    pcolor(lead,varnum,varplot)
    set(gca,'color','w')
    caxis([0.55 1])
    cmocean('amp',9)
    
    subplot(4,ni,ii+4),hold on
    varplot = nan(nv+1,nl+1);
    varplot(1:end-1,1:end-1) = squeeze(sbias(ii,:,:));
    pcolor(lead,varnum,varplot)
    caxis([-2 2])
    cmocean('balance',10)
        
    subplot(4,ni,ii+6),hold on
    varplot = nan(nv+1,nl+1);
    varplot(1:end-1,1:end-1) = squeeze(srmsd(ii,:,:));
    pcolor(lead,varnum,varplot)
    caxis([0 2])
    cmocean('amp',10)
end

% Format plots
for ii = 1:ni*4
    subplot(4,ni,ii)
    axis([0 nl 0 nv])
    set(gca,'fontsize',12,'tickdir','out')
    set(gca,'xtick',0:12,'ytick',.5:nv-.5)
    if ii>3*ni
        set(gca,'xticklabel',0:12)
        xlabel('Lead time (months)')
    else
        set(gca,'xticklabel',[])
    end
    if rem(ii,ni)==1
        set(gca,'yticklabel',vars)
    else
        set(gca,'yticklabel',[])
    end
    set(gca,'position',[xpos(ii) ypos(ii) width height])
    if rem(ii,ni)==0
        cb = colorbar;
        set(cb,'fontsize',12,'tickdir','out')
        set(cb,'position',[xpos(ii)+width+.02 ypos(ii) .02 height])
        if ii==ni
            ylabel(cb,'ACC')
        elseif ii==2*ni
            ylabel(cb,'Forecast Accuracy')
        elseif ii==3*ni
            ylabel(cb,'Standardized Bias')
        elseif ii==4*ni
            ylabel(cb,'Standardized RMSD')
        end
    end
    if ii==1
        title('Jan initialization')
    elseif ii==2
        title('Jul initialization')
    end
end

% ==========================
function fig3
% ==========================
% Skill maps for all variables

% Load input file
f_in = 'fig3.mat';
load(f_in)

% Plotting info
xpos = [.03 .31 .59];
ypos = [.69 .37 .05];
width = .26;
height = .26;
cbar_pos = [xpos(end)+width+.02 ypos(end) .02 ypos(1)+height-ypos(end)];

% Open figure
figure,set(gcf,'color','w')
set(gcf,'position',[50 50 450 420])

% Loop through variables and plot
for iv = 1:length(vars)
    subplot(3,3,iv),hold on
    pcolor(lon,lat,squeeze(var_plot(iv,:,:)));
    title(vars{iv})
end

% Format plot
% Modified colorbar with white in middle
ctmp = cmocean('balance',20);
ctmp(10:11,:) = 1;
colormap(ctmp)
kk = 1;
for ii = 1:3
    for jj = 1:3
        subplot(3,3,kk)
        hold on
        shading interp
        set(gca,'color',[.8 .8 .8])
        set(gca,'fontsize',12,'xticklabels',[],'yticklabels',[])
        caxis([-1 1])
        box on
        axis([-128 -116 31 47])
        plot(coastline(:,1),coastline(:,2),'k.')
        set(gca,'position',[xpos(jj) ypos(ii) width height])
        if kk==6
            cb = colorbar;
            set(cb,'tickdir','out','fontsize',12,'position',cbar_pos)
            ylabel(cb,'ACC')
        end
        kk = kk+1;
    end
end

% ==========================
function fig4
% ==========================
% Skill maps for SST

% Load input file
f_in = 'fig4.mat';
load(f_in)

% Plotting info
figpos = [100 100 450 400];
xpos = [.06 .48 .06 .48];
ypos = [.52 .52 .05 .05];
width = .38;
height = .43;
cbar_pos = [xpos(end)+width+.02 ypos(end) .02 ypos(1)+height-ypos(end)];

% Open figure
figure,set(gcf,'color','w')
set(gcf,'position',figpos)

% Plot ACC
subplot(221),hold on
pcolor(lon,lat,acc_jan(:,:,1));
ind = find(sig_zero_jan(:,:,1)==0 & rem(lon,.2)==0 & rem(lat,.2)==0);
plot(lon(ind),lat(ind),'k.','markersize',4)
title('1.5-month lead')
ylabel('January initialization','fontweight','bold')

subplot(222),hold on
pcolor(lon,lat,acc_jan(:,:,2));
ind = find(sig_zero_jan(:,:,2)==0 & rem(lon,.2)==0 & rem(lat,.2)==0);
plot(lon(ind),lat(ind),'k.','markersize',4)
title('6.5-month lead')

subplot(223),hold on
pcolor(lon,lat,acc_jul(:,:,1));
ind = find(sig_zero_jul(:,:,1)==0 & rem(lon,.2)==0 & rem(lat,.2)==0);
plot(lon(ind),lat(ind),'k.','markersize',4)
ylabel('July initialization','fontweight','bold')

subplot(224),hold on
pcolor(lon,lat,acc_jan(:,:,2));
ind = find(sig_zero_jan(:,:,2)==0 & rem(lon,.2)==0 & rem(lat,.2)==0);
plot(lon(ind),lat(ind),'k.','markersize',4)

% Format plot
for ii = 1:4
    subplot(2,2,ii)
    shading interp
    set(gca,'color',[.8 .8 .8])
    set(gca,'fontsize',12,'xticklabels',[],'yticklabels',[])
    caxis([-1 1])
    box on
    axis([-130 -116 30 48])
    plot(coastline(:,1),coastline(:,2),'k.')
    set(gca,'position',[xpos(ii) ypos(ii) width height])

    % Modified colorbar with white in middle
    ctmp = cmocean('balance',20);
    ctmp(10:11,:) = 1;
    colormap(ctmp)
    if ii==2
        cb = colorbar;
        set(cb,'tickdir','out','fontsize',12,'position',cbar_pos)
        ylabel(cb,'ACC')
    end
end

% ==========================
function fig5
% ==========================
% Skill maps for SSH

% Load input file
f_in = 'fig5.mat';
load(f_in)

% Plotting info
figpos = [100 100 450 400];
xpos = [.06 .48 .06 .48];
ypos = [.52 .52 .05 .05];
width = .38;
height = .43;
cbar_pos = [xpos(end)+width+.02 ypos(end) .02 ypos(1)+height-ypos(end)];

% Open figure
figure,set(gcf,'color','w')
set(gcf,'position',figpos)

% Plot ACC
subplot(221),hold on
pcolor(lon,lat,acc_jan(:,:,1));
ind = find(sig_zero_jan(:,:,1)==0 & rem(lon,.2)==0 & rem(lat,.2)==0);
plot(lon(ind),lat(ind),'k.','markersize',4)
title('1.5-month lead')
ylabel('January initialization','fontweight','bold')

subplot(222),hold on
pcolor(lon,lat,acc_jan(:,:,2));
ind = find(sig_zero_jan(:,:,2)==0 & rem(lon,.2)==0 & rem(lat,.2)==0);
plot(lon(ind),lat(ind),'k.','markersize',4)
title('6.5-month lead')

subplot(223),hold on
pcolor(lon,lat,acc_jul(:,:,1));
ind = find(sig_zero_jul(:,:,1)==0 & rem(lon,.2)==0 & rem(lat,.2)==0);
plot(lon(ind),lat(ind),'k.','markersize',4)
ylabel('July initialization','fontweight','bold')

subplot(224),hold on
pcolor(lon,lat,acc_jul(:,:,2));
ind = find(sig_zero_jul(:,:,2)==0 & rem(lon,.2)==0 & rem(lat,.2)==0);
plot(lon(ind),lat(ind),'k.','markersize',4)

% Format plot
for ii = 1:4
    subplot(2,2,ii)
    shading interp
    set(gca,'color',[.8 .8 .8])
    set(gca,'fontsize',12,'xticklabels',[],'yticklabels',[])
    caxis([-1 1])
    box on
    axis([-130 -116 30 48])
    plot(coastline(:,1),coastline(:,2),'k.')
    set(gca,'position',[xpos(ii) ypos(ii) width height])

    % Modified colorbar with white in middle
    ctmp = cmocean('balance',20);
    ctmp(10:11,:) = 1;
    colormap(ctmp)
    if ii==2
        cb = colorbar;
        set(cb,'tickdir','out','fontsize',12,'position',cbar_pos)
        ylabel(cb,'ACC')
    end
end

% ==========================
function fig6(panel)
% ==========================
% Skill maps for BT

% Load input file
% Each panel of figure is plotted separately
% Replace fig6a with fig6b, fig6c, fig6d for other panels
f_in = sprintf('fig6%s.mat',panel);
load(f_in)

% Plotting info
% load ~/Documents/Data/CCSRA/ccsra_info h coastline
% max_depth = 1000; % Max depth to plot
% %axlims = [-126 -116 30 48];
axlims = [-127.5 -121.5 40 48;...
          -126 -120 35 40;...
          -121.5 -115.5 30 35];
xpos = 0.1;
width = 0.6;
tot_height = 0.95;
height = (axlims(:,4)-axlims(:,3))/18*tot_height;
yspace = 0.01;
ypos(3) = .02;
ypos(2) = ypos(3)+height(3)+yspace;
ypos(1) = ypos(2)+height(2)+yspace;
    
% Modified colorbar with white in middle
ctmp = cmocean('balance',20);
ctmp(10:11,:) = 1;

% % Input directory
% dirin = '~/Documents/Data/CanCM4/CanCM4_ROMS/forecasts';
% 
% % Load skill metrics
% fname = sprintf('%s/bt_skill_metrics_map',dirin);
% load(fname,'skill')
% 
% % Plot ACC
% if init_mon==1
%     forc_plot = squeeze(skill.bc.acc(1,:,:,end,lead_plot));
%     sig_zero_plot = squeeze(skill.bc.acc_sig_zero(1,:,:,end,lead_plot));
% elseif init_mon==7
%     forc_plot = squeeze(skill.bc.acc(2,:,:,end,lead_plot));
%     sig_zero_plot = squeeze(skill.bc.acc_sig_zero(2,:,:,end,lead_plot));
% else
%     disp('init_mon must be 1 or 7')
% end
% 
% % Mask out deeper areas
% forc_plot(h>max_depth) = 0;
% sig_zero_plot(h>max_depth) = nan;

% Open figure
figure,set(gcf,'color','w'),hold on
set(gcf,'position',[100 100 200 350])

% Plot
for ii = [3 2 1]
    subplot(3,1,ii),hold on
    pcolor(lon-.05,lat-.05,squeeze(acc));
    shading flat
    ind = find(squeeze(sig_zero)==0);
    plot(lon(ind),lat(ind),'.','markersize',4,'color',[.2 .2 .2])
    caxis([-1 1])
    colormap(ctmp)
    plot(coastline(:,1),coastline(:,2),'k.','markersize',2)
    axis(axlims(ii,:))
    set(gca,'xticklabel',[],'tickdir','out')
    set(gca,'color',[.8 .8 .8])
    if ii==3
        set(gca,'ytick',axlims(ii,3):axlims(ii,4))
    else
        set(gca,'ytick',axlims(ii,3)+1:axlims(ii,4))
    end
    if ii==2
        cb = colorbar;
        set(cb,'fontsize',12,'tickdir','out')
        ylabel(cb,'ACC')
        set(cb,'position',[xpos+width+.02 ypos(3) .04 tot_height+2*yspace])
    end
    set(gca,'position',[xpos ypos(ii) width height(ii)])
    box on
end
 
% ==========================
function fig7
% ==========================
% Evaluate SST skill with shore stations

% Load input file
f_in = 'fig7.mat';
load(f_in)

% Plotting info
cplot = [.7 0 .7;1 .6 0;0 .7 0]; % Colors for plotting
xpos = [.05 .58 .8];
ypos = [.7:-.2:.1];
width = [.45 .18];
height(1) = .17;
height(2) = ypos(1)-ypos(end)+height(1); 
ylim = [4 4 4 4];

% Plot time series
figure,set(gcf,'color','w','position',[100 100 800 600])
ns = length(station_lat);
ny = length(year_forc);
for is = 1:ns
    subplot(ns,3,3*is-2),hold on
    h1(2) = plot(decyear_obs,ssta_obs(is,:),'color','k','linewidth',2);
    
    for iy = 1:ny
        % Make monthly time coordinate to plot each year
        decyear_jan = year_forc(iy)+1/24:1/12:year_forc(iy)+23/24;
        
        % Plot jan initialized forecasts
        forc_plot = squeeze(ssta_forc(is,:,iy));
        h1(1) = plot(decyear_jan,forc_plot,'color',cplot(3,:),'linewidth',2);
        
        % Plot forecast starts
        plot(decyear_jan(1),forc_plot(1),'o','color',cplot(3,:),'markerfacecolor',cplot(3,:),'markersize',7)
    end
    
    % Format plot
    set(gca,'tickdir','out','fontsize',12)
    box on,grid on
    if is<ns
        set(gca,'xticklabel',[])
    else
        xlabel('Year')
    end
    if is==round(ns/2)
        ylabel(['SST anomaly (' char(176) 'C)'])
    end
    if is==1
        hl1 = legend(h1,'Jan Forecast','Observed');
        hl1_width = hl1.Position(3);
        xspace = width(1) - hl1_width; % Used to center legend
        set(hl1,'fontsize',12,'orientation','horizontal','position',[xpos(1)+xspace/2 ypos(1)+height(1)+.02 hl1_width .03])
    end
    axis([1982 2011 -ylim(is) ylim(is)])
    text(0.01,0.1,station_name{is},'units','normalized','fontsize',12)
end

% Plot correlations
subplot(ns,3,2),hold on
for il = 1:2
    h2(il) = plot(squeeze(acc_forc_jan(:,il)),station_lat,'color',cplot(il,:),'linewidth',2);
    plot(acc_pers_jan(:,il),station_lat,'--','color',cplot(il,:),'linewidth',2);
    plot(acc_pers_jan(:,il),station_lat,'o','color',cplot(il,:),'markerfacecolor','w','linewidth',2)
    plot(squeeze(acc_forc_jan(:,il)),station_lat,'o','color',cplot(il,:),'markerfacecolor',cplot(il,:),'linewidth',2)
end
set(gca,'tickdir','out','fontsize',12)
set(gca,'xtick',-.5:.5:1,'xticklabel',-.5:.5:1,'ytick',31:2:47,'yticklabel',31:2:47)
axis([-.5 1 32 38])
grid on,box on
ylabel('Latitude')
xlabel('ACC')
title('Jan Initialization')

subplot(ns,3,3),hold on
for il = 1:2
    plot(squeeze(acc_forc_jul(:,il)),station_lat,'color',cplot(il,:),'linewidth',2);
    plot(acc_pers_jul(:,il),station_lat,'--','color',cplot(il,:),'linewidth',2);
    plot(acc_pers_jul(:,il),station_lat,'o','color',cplot(il,:),'markerfacecolor','w','linewidth',2);
    plot(squeeze(acc_forc_jul(:,il)),station_lat,'o','color',cplot(il,:),'markerfacecolor',cplot(il,:),'linewidth',2);
end
set(gca,'tickdir','out','fontsize',12)
set(gca,'xtick',-.5:.5:1,'xticklabel',-.5:.5:1,'ytick',31:2:47,'yticklabel',[])
axis([-.5 1 32 38])
grid on,box on
xlabel('ACC')
title('Jul Initialization')

% Legend
h2(il+1) = plot(1e6,1e6,'--','color',[.3 .3 .3],'linewidth',2); % Dummy plot
hl2 = legend(h2,'1.5-month lead','6.5-month lead','Persistence');
set(hl2,'fontsize',12,'location','southwest')

% Set panel positions
subplot(ns,3,3)
set(gca,'position',[xpos(3) ypos(end) width(2) height(2)])
subplot(ns,3,2)
set(gca,'position',[xpos(2) ypos(end) width(2) height(2)])
for is = ns:-1:1
    subplot(ns,3,3*is-2)
    set(gca,'position',[xpos(1) ypos(is) width(1) height(1)])
end

% ==========================
function fig8
% ==========================
% Evaluate SSH skill with shore stations

% Load input file
f_in = 'fig8.mat';
load(f_in)

% Plotting info
cplot = [.7 0 .7;1 .6 0;0 .7 0]; % Colors for plotting
xpos = [.07 .59 .81];
ypos = [.87:-.1:.07];
width = [.45 .18];
height(1) = .085;
height(2) = ypos(1)-ypos(end)+height(1); 
ylim = [.4 .4 .4 .3 .3 .2 .2 .2 .2];

% Plot time series
figure,set(gcf,'color','w','position',[100 100 800 600])
ns = length(station_lat);
ny = length(year_forc);
for is = 1:ns
    subplot(ns,3,3*is-2),hold on
    h1(2) = plot(decyear_obs,ssha_obs(is,:),'color','k','linewidth',2);
    
    for iy = 1:ny
        % Make monthly time coordinate to plot each year
        decyear_jan = year_forc(iy)+1/24:1/12:year_forc(iy)+23/24;
        
        % Plot jan initialized forecasts
        forc_plot = squeeze(ssha_forc(is,:,iy));
        h1(1) = plot(decyear_jan,forc_plot,'color',cplot(3,:),'linewidth',2);
        
        % Plot forecast starts
        plot(decyear_jan(1),forc_plot(1),'o','color',cplot(3,:),'markerfacecolor',cplot(3,:),'markersize',7)
    end
    
    % Format plot
    set(gca,'tickdir','out','fontsize',12)
    box on,grid on
    if is<ns
        set(gca,'xticklabel',[])
    else
        xlabel('Year')
    end
    if is==round(ns/2)
        ylabel(['SSH anomaly (m)'])
    end
    if is==1
        %hl1 = legend(h1,'Jan Initialization','Jul Initialization','Observed');
        hl1 = legend(h1,'Jan Initialization','Observed');
        hl1_width = hl1.Position(3);
        xspace = width(1) - hl1_width; % Used to center legend
        set(hl1,'fontsize',12,'orientation','horizontal','position',[xpos(1)+xspace/2 ypos(1)+height(1)+.01 hl1_width .03])
    end
    axis([1982 2011 -ylim(is) ylim(is)])
    text(0.01,0.1,station_name{is},'units','normalized','fontsize',12)
end

% Plot correlations
subplot(ns,3,2),hold on
for il = 1:2
    h2(il) = plot(squeeze(acc_forc_jan(:,il)),station_lat,'color',cplot(il,:),'linewidth',2);
    plot(acc_pers_jan(:,il),station_lat,'--','color',cplot(il,:),'linewidth',2);
    plot(acc_pers_jan(:,il),station_lat,'o','color',cplot(il,:),'markerfacecolor','w','linewidth',2)
    plot(squeeze(acc_forc_jan(:,il)),station_lat,'o','color',cplot(il,:),'markerfacecolor',cplot(il,:),'linewidth',2)
end
set(gca,'tickdir','out','fontsize',12)
set(gca,'xtick',0:.5:1,'xticklabel',0:.5:1,'ytick',31:2:47,'yticklabel',31:2:47)
axis([0 1 31 48])
grid on,box on
ylabel('Latitude')
xlabel('ACC')
title('Jan Initialization')

subplot(ns,3,3),hold on
for il = 1:2
    plot(squeeze(acc_forc_jul(:,il)),station_lat,'color',cplot(il,:),'linewidth',2);
    plot(acc_pers_jul(:,il),station_lat,'--','color',cplot(il,:),'linewidth',2);
    plot(acc_pers_jul(:,il),station_lat,'o','color',cplot(il,:),'markerfacecolor','w','linewidth',2);
    plot(squeeze(acc_forc_jul(:,il)),station_lat,'o','color',cplot(il,:),'markerfacecolor',cplot(il,:),'linewidth',2);
end
set(gca,'tickdir','out','fontsize',12)
set(gca,'xtick',0:.5:1,'xticklabel',0:.5:1,'ytick',31:2:47,'yticklabel',[])
axis([0 1 31 48])
grid on,box on
xlabel('ACC')
title('Jul Initialization')

% Legend
h2(il+1) = plot(1e6,1e6,'--','color',[.3 .3 .3],'linewidth',2); % Dummy plot
hl2 = legend(h2,'1.5-month lead','6.5-month lead','Persistence');
set(hl2,'fontsize',12,'location','southwest')

% Set panel positions
subplot(ns,3,3)
set(gca,'position',[xpos(3) ypos(end) width(2) height(2)])
subplot(ns,3,2)
set(gca,'position',[xpos(2) ypos(end) width(2) height(2)])
for is = ns:-1:1
    subplot(ns,3,3*is-2)
    set(gca,'position',[xpos(1) ypos(is) width(1) height(1)])
end

% ==========================
function fig9
% ==========================
% Compare forecast ensemble spread to observed variability

% Load input file
f_in = 'fig9.mat';
load(f_in)

% Plotting info
clim = [0.5 1.5;0 .1;0 1.2;0 3e-3;0.6 1.1;.03 .08;0 1;.2e-3 2e-3];
ncol = [10 10 12 12 10 10 10 12];
xpos = [.05 .29 .53 .77];
ypos = [.53 .05];
width = .17;
height = .43;
titles = {['SST (' char(176) 'C)'] 'SSH (m)' ['BT (' char(176) 'C)'] 'BF (s^{-1})'};
xpos = [xpos xpos];
ypos = [ypos(1)*ones(1,4) ypos(2)*ones(1,4)];

% Open figure
figure,set(gcf,'color','w','position',[50 50 800 500])

% Plot
nv = length(vars);
for iv = 1:nv
    subplot(2,nv,iv),hold on
    pcolor(lon-.05,lat-.05,squeeze(spread(iv,:,:)))
    
    subplot(2,nv,iv+nv),hold on
    pcolor(lon-.05,lat-.05,squeeze(std_obs(iv,:,:)))
end

% Format panels
for ii = 1:2*nv
    subplot(2,nv,ii)
    if ismember(ii,[1 2 5 6])
        axis([-132 -116 30.5 47.5])
    else
        axis([-125 -116 30.5 47.5])
    end
    set(gca,'tickdir','out','fontsize',10,'color',[.9 .9 .9])
    if ii<5
        set(gca,'xticklabels',[])
    end
    if ~ismember(ii,[1 5])
        set(gca,'yticklabels',[])
    end
    caxis(clim(ii,:))
    cmocean('thermal',ncol(ii))
    plot(coastline(:,1),coastline(:,2),'k.','markersize',2)
    box on
    shading flat
    set(gca,'position',[xpos(ii) ypos(ii) width height])
    if ii<=length(titles)
        title(titles{ii},'fontsize',12)
    end
    if ii==1
        ylabel('Model spread','fontweight','bold','fontsize',12)
    elseif ii==nv+1
        ylabel('Observed variability','fontweight','bold','fontsize',12)
    end
    cb = colorbar;
    if ncol(ii)==10
        set(cb,'ytick',linspace(clim(ii,1),clim(ii,2),6))
    elseif ncol(ii)==12
        set(cb,'ytick',linspace(clim(ii,1),clim(ii,2),7))
    end
    set(cb,'tickdir','out','fontsize',10)
    set(cb,'position',[xpos(ii)+width+.005 ypos(ii) .015 height])
end

% ==========================
function fig10
% ==========================
% Skill summary for individual ensemble members

% Load input file
f_in = 'fig10.mat';
load(f_in)

% Bar plot of mean skill for each variable/member averaged across all
% leads/initializations
figure,set(gcf,'color','w','position',[50 50 450 300])
hold on
x = categorical(vars);
x = reordercats(x,vars);
h = bar(x,acc_mean);
cmap = cmocean('thermal',7);
h(1).FaceColor = cmap(2,:);
h(2).FaceColor = cmap(4,:);
h(3).FaceColor = cmap(6,:);
h(4).FaceColor = [.2 .2 .2];
h(1).BarWidth = 1;
set(gca,'fontsize',12,'tickdir','out')
ylabel('ACC')
hl = legend('Member 2','Member 8','Member 10','Ensemble mean');
set(hl,'fontsize',12)
box on

% ==========================
function fig11
% ==========================
% Skill summary for all variables and metrics
% Compare interpolated and bias corrected forcing for each initialization
% month

% Load input file
f_in = 'fig11.mat';
load(f_in)
 
% Plotting info
xpos = [.08 .5 .08 .5 .08 .5 .08 .5]; % For 2 columns
ypos = [.75 .75 .53 .53 .31 .31 .09 .09];
width = .36;
height = .19;

% Open figure
figure,set(gcf,'color','w','position',[100 100 450 750])

% Plot
[~,nv,nl] = size(acc_bc);

% Lead and variable matrices for plotting
leadmat = repmat(0:12,nv+1,1);
varmat = repmat([0:nv]',1,nl+1);

subplot(4,2,1),hold on
varplot = nan(nv+1,nl+1); % Pad for plotting with pcolor
varplot(1:end-1,1:end-1) = squeeze(acc_bc(1,:,:)-acc_int(1,:,:));
pcolor(leadmat,varmat,varplot)
caxis([-1 1])

subplot(4,2,2),hold on
varplot = nan(nv+1,nl+1);
varplot(1:end-1,1:end-1) = squeeze(acc_bc(2,:,:)-acc_int(2,:,:));
pcolor(leadmat,varmat,varplot)
caxis([-1 1])

subplot(4,2,3),hold on
varplot = nan(nv+1,nl+1);
varplot(1:end-1,1:end-1) = squeeze(fa_bc(1,:,:)-fa_int(1,:,:));
pcolor(leadmat,varmat,varplot)
caxis([-1 1])

subplot(4,2,4),hold on
varplot = nan(nv+1,nl+1);
varplot(1:end-1,1:end-1) = squeeze(fa_bc(2,:,:)-fa_int(2,:,:));
pcolor(leadmat,varmat,varplot)
caxis([-1 1])

subplot(4,2,5),hold on
varplot = nan(nv+1,nl+1);
varplot(1:end-1,1:end-1) = squeeze(abs(sbias_bc(1,:,:))-abs(sbias_int(1,:,:)));
pcolor(leadmat,varmat,varplot)
caxis([-2 2])

subplot(4,2,6),hold on
varplot = nan(nv+1,nl+1);
varplot(1:end-1,1:end-1) = squeeze(abs(sbias_bc(2,:,:))-abs(sbias_int(2,:,:)));
pcolor(leadmat,varmat,varplot)
caxis([-2 2])

subplot(4,2,7),hold on
varplot = nan(nv+1,nl+1); 
varplot(1:end-1,1:end-1) = squeeze(srmsd_bc(1,:,:)-srmsd_int(1,:,:));
pcolor(leadmat,varmat,varplot)
caxis([-2 2])

subplot(4,2,8),hold on
varplot = nan(nv+1,nl+1);
varplot(1:end-1,1:end-1) = squeeze(srmsd_bc(2,:,:)-srmsd_int(2,:,:));
pcolor(leadmat,varmat,varplot)
caxis([-2 2])


% Format plots
for ii = 1:8
    subplot(4,2,ii)
    axis([0 nl 0 nv])
    set(gca,'fontsize',12,'tickdir','out')
    set(gca,'xtick',0:12,'ytick',.5:nv-.5)
    if ii>6
        set(gca,'xticklabel',0:12)
        xlabel('Lead time (months)')
    else
        set(gca,'xticklabel',[])
    end
    if rem(ii,2)==1
        set(gca,'yticklabel',vars)
    else
        set(gca,'yticklabel',[])
    end
    set(gca,'position',[xpos(ii) ypos(ii) width height])
    cmocean('balance',20)
    if rem(ii,2)==0
        cb = colorbar;
        set(cb,'fontsize',12,'tickdir','out')
        set(cb,'position',[xpos(ii)+width+.02 ypos(ii) .02 height])
        if ii==2
            ylabel(cb,'\Delta ACC')
        elseif ii==4
            ylabel(cb,'\Delta Forecast Accuracy')
        elseif ii==6
            ylabel(cb,'\Delta Standardized Bias')
        elseif ii==8
            ylabel(cb,'\Delta Standardized RMSD')
        end
    end
    if ii==1
        title('Jan Initialization')
    elseif ii==2
        title('Jul Initialization')
    end
end

% ==========================
function fig12
% ==========================
% Skill summary for LMR applications
%
% Application summary
% #   NAME          VARS    REGIONS     MONTHS
% --  ----          ----    -------     ------
% 1   ecocast       all     all         Aug-Jan
% 2   total         sst     south       Dec-Jul
% 3   WhaleWatch    all     south       May-Nov
% 4   HCI           sst     central     Mar-Nov
% 5   sardine       sst,ssh south       Mar-May
% 6   albacore      all     north       Jun-Nov
% 7   AEI           all     central     Apr-Jul

% Load input file
f_in = 'fig12.mat';
load(f_in)
    
% Color for plotting
cplot = [0 .4 0;0 .6 0]; % Green

% Box plot
figure,hold on
set(gcf,'position',[50 50 500 450],'color','w')
plot([0 0],[0 21],'k:')
boxplot(acc_compiled','orientation','horizontal','widths',.6,...
    'color',[cplot(2,:);cplot(1,:);cplot(1,:)],'symbol','')

% Plot mean
plot(acc_mean(1:3:19),1:3:19,'o','color',cplot(2,:),'markerfacecolor',cplot(2,:))
plot(acc_mean(2:3:20),2:3:20,'o','color',cplot(1,:),'markerfacecolor',cplot(1,:))

% Format boxplot
h = findobj(gca,'tag','Box');
set(h(1:3:19),'linewidth',1.5)
h = findobj(gca,'tag','Lower Whisker');
set(h,'linestyle','-')
set(h(1:3:19),'linewidth',1.5,'linestyle','-')
h = findobj(gca,'tag','Upper Whisker');
set(h,'linestyle','-')
set(h(1:3:19),'linewidth',1.5)
h = findobj(gca,'tag','Lower Adjacent Value');
set(h(1:3:19),'linewidth',1.5)
h = findobj(gca,'tag','Upper Adjacent Value');
set(h(1:3:19),'linewidth',1.5)
h = findobj(gca,'tag','Median');
set(h(1:3:19),'linewidth',1.5)

% Format plot
set(gca,'ytick',2:3:20,'yticklabel',apps)
set(gca,'tickdir','out','fontsize',12,'xtick',-.4:.2:1)
xlabel('ACC')
ylabel('Application')
axis([-.5 1 0 21])

% Legend
clear h
h(1) = fill([-10 -11 -11 -10],[0 0 1 1],'w','linewidth',1.5,'edgecolor',cplot(1,:));
h(2) = fill([-10 -11 -11 -10],[0 0 1 1],'w','linewidth',.5,'edgecolor',cplot(2,:));
hl = legend(h,'Full model','Reduced model');
set(hl,'fontsize',12,'location','northwest')

% ==========================
function fig13
% ==========================
% Compare model forecast skill to persistence

% Load input file
f_in = 'fig13.mat';
load(f_in)

% Plotting info
titles = {'North' 'Central' 'South'};
xpos = .13;
ypos = [.71 .39 .07];
width = .83;
height = .26;

% % Calculate skill
% skill_sst = roms_skill_metrics('sst',dist_bnds);
% skill_ssh = roms_skill_metrics('ssh',dist_bnds);
% skill_svstr = roms_skill_metrics('svstr',dist_bnds);
% skill_sustr = roms_skill_metrics('sustr',dist_bnds);
% 
% % Choose initialization month
% if init_mon==1
%     month = 'January';
% elseif init_mon==7
%     month = 'July';
%     init_mon = 2;
% else
%     disp('init_mon must be 1 or 7')
%     return
% end

% Open figure
figure,set(gcf,'color','w','position',[50 50 400 600])


% Loop through regions and plot skill above persistence
jj = 1;
for ii = 3:-1:1
%     % Load skill
%     acc_sst = squeeze(skill_sst.bc.acc(init_mon,ii,end,:));
%     acc_sst_pers = squeeze(skill_sst.pers.acc(init_mon,ii,:));
%     acc_ssh = squeeze(skill_ssh.bc.acc(init_mon,ii,end,:));
%     acc_ssh_pers = squeeze(skill_ssh.pers.acc(init_mon,ii,:));
%     acc_svstr = squeeze(skill_svstr.bc.acc(init_mon,ii,end,:));
%     acc_svstr_pers = squeeze(skill_svstr.pers.acc(init_mon,ii,:));
%     acc_sustr = squeeze(skill_sustr.bc.acc(init_mon,ii,end,:));
%     acc_sustr_pers = squeeze(skill_sustr.pers.acc(init_mon,ii,:));
%     
%     % Calculate correlations
%     r_ssh = corr(acc_sst-acc_sst_pers,acc_ssh-acc_ssh_pers);
%     r_wind = corr(acc_sst(2:end)-acc_sst_pers(2:end),(acc_svstr(1:end-1)-acc_svstr_pers(1:end-1)+acc_sustr(1:end-1)-acc_sustr_pers(1:end-1))/2);

    % Plot
    subplot(3,1,jj),hold on
    plot([0 12],[0 0],'k--','linewidth',.5)
    h(1) = plot(.5:11.5,acc_diff_sst(ii,:),'k','linewidth',2);
    h(2) = plot(1.5:12.5,acc_diff_wind(ii,:),'b','linewidth',2);
    h(3) = plot(.5:11.5,acc_diff_ssh(ii,:),'r','linewidth',2);
    text(.45,.08,sprintf('r = %.2f',r_wind(ii)),'color','b','backgroundcolor','w','fontsize',12,'fontweight','bold','units','normalized')
    text(.62,.08,sprintf('r = %.2f',r_ssh(ii)),'color','r','backgroundcolor','w','fontsize',12,'fontweight','bold','units','normalized')
    
    % Format plot
    axis([0 12 -.5 1])
    set(gca,'tickdir','out','xtick',0:12,'fontsize',12)
    grid on,box on
    if jj==1
        hl = legend(h,'sst','wind (1-month lead)','ssh');
        set(hl,'location','northwest','fontsize',12)
    elseif jj==2
        ylabel('\DeltaACC')
    end
    if jj<3
        set(gca,'xticklabel',[])
    else
        xlabel('Lead time (months)')
    end
    title(titles{jj})
    set(gca,'position',[xpos ypos(jj) width height])
    jj = jj+1;
end

% ==========================
function fig14
% ==========================
% SSH and SST skill maps with verification from sensitivity runs

% Load input file
f_in = 'fig14.mat';
load(f_in)

% Plotting info
xpos = [.04 .26 .48 .7];
ypos = [.55 .1];
width = .2;
height = .4;
varlabels = {'SSH' 'SST'};

% Runs
run_title = {'CCSRA' 'Hindcast' 'Wind' 'Ocean'};
nr = length(run_title);

% Open figure
figure,set(gcf,'color','w','position',[50 50 900 500])

% Make colormap with 20 bins
cmap = cmocean('balance',20);
colormap(cmap)

% Plot
load ~/Documents/Data/CCSRA/ccsra_info coastline
for ii = 1:nr
    subplot(2,nr,ii),hold on
%     cmd = sprintf('tmp = acc_ssh.bc.%s;',run_name{ii});
%     eval(cmd)
    pcolor(lon,lat,acc_ssh(:,:,ii))
    
    subplot(2,nr,ii+nr),hold on
%     cmd = sprintf('tmp = acc_sst.bc.%s;',run_name{ii});
%     eval(cmd)
    pcolor(lon,lat,acc_sst(:,:,ii))
end
    
% Format plot
ii = 1;
for iv = 1:2
    for ir = 1:nr
        subplot(2,nr,ii)
        shading interp
        set(gca,'color',[.8 .8 .8])
        set(gca,'tickdir','out','fontsize',12,'xticklabel',[],'yticklabel',[])
        caxis([-1 1])
        axis([-130 -116 30 48])
        box on
        plot(coastline(:,1),coastline(:,2),'k.')
        if iv==1
            title(run_title{ir})
        end
        if ir==1
            ylabel(varlabels{iv},'fontweight','bold')
        end
        set(gca,'position',[xpos(ir) ypos(iv) width height])
        if ir==nr & iv==2
            cb = colorbar;
            set(cb,'tickdir','out','fontsize',12)
            set(cb,'position',[xpos(end)+width+.01 ypos(2) .02 ypos(1)+height-ypos(2)])
            ylabel(cb,'ACC')
        end
        ii = ii+1;
    end
end
