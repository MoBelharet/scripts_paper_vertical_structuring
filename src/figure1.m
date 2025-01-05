
project = 'MALASPINA';


load('lat_sa_m.mat')
load('lon_sa_m.mat')

d = importdata('cluster_list_selected.txt');
clstr_list = d.data;

lat = lat_sa_m(clstr_list);
lon = lon_sa_m(clstr_list);


load("id_.mat")
id_(37:38) = [];
load("id_over_est.mat")
id_over_est = [84:93, 95];
load("id_under_est.mat")
%% Figure 1-a

col_id = {[75 0 130]/225,[0 1 0],[1 0 0]};

fig = figure('Renderer', 'painters', 'Position', [50 5 1600 1200]);

ax1 = subplot(2,16,1:16);

worldmap('World');

axesm ('pcarree', 'Frame', 'on', 'Grid', 'on');
setm(ax1, 'Origin', [],'MapLatLimit',[-50 50],'MeridianLabel','on','ParallelLabel','on')
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax1, land, 'FaceColor', [0.9 0.9 0.9])

set(findall(ax1,'Tag','MLabel'),'visible','on','fontsize',10) % x-axis
set(findall(ax1,'Tag','PLabel'),'visible','on','fontsize',10) % y_axis

scatterm(ax1,lat(id_),lon(id_),50,col_id{1},'o','filled','MarkerEdgeColor','k')
scatterm(ax1,lat(id_over_est),lon(id_over_est),50,col_id{2},'o','filled', 'MarkerEdgeColor','k')
scatterm(ax1,lat(id_under_est),lon(id_under_est),50,col_id{3},'o','filled', 'MarkerEdgeColor','k')

%scatterm(ax1,lat(id_(37:38)),lon(id_(37:38)),50,'w','o','filled','MarkerEdgeColor','w','markerfacealpha',1)
 
scatterm(ax1,-26, -165,100,col_id{1},'o','filled','MarkerEdgeColor','k')
scatterm(ax1,-34, -165,100,col_id{2},'o','filled','MarkerEdgeColor','k')
scatterm(ax1,-42, -165,100,col_id{3},'o','filled','MarkerEdgeColor','k')


textm(-26, -160, 'GP stations','fontsize',13,'color',col_id{1})
textm(-34, -160, 'SP stations','fontsize',13,'color',col_id{2})
textm(-42, -160, 'DP stations','fontsize',13,'color',col_id{3})

%%
ax2 = subplot(2,16,19:23);
ax3 = subplot(2,16,26:30);

%%
id_depth = 1:100;
    

   
load(['profile_total_' project '_' prj ])

load('depth')
id = find(depth<=1000);
depth = depth(id);
if(size(depth,1)>1)
    depth = depth';
end
depth = -depth(id_depth);

d = importdata('cluster_list_selected.txt');
clstr = d.data;

load('day_sa_m')
load('night_sa_m')

day_sa = day_sa_m(id,clstr);
night_sa = night_sa_m(id,clstr);

ind_d = find(isnan(day_sa));
ind_n = find(isnan(night_sa));


load('oxy_obs_m')
oxy_ = nanmin(oxy_obs_m(1:100,clstr));


profile_total_d = squeeze(nanmean(profile_total(1,:,:,:),4));
profile_total_n = squeeze(nanmean(profile_total(2,:,:,:),4));

DEPTH = repmat(depth,2,1,length(clstr));
 
profile_total_d(ind_d) = nan;
profile_total_n(ind_n) = nan;

profile_total_d = profile_total_d(id_depth,:)./(nansum(profile_total_d(id_depth,:)));
profile_total_n = profile_total_n(id_depth,:)./(nansum(profile_total_n(id_depth,:)));

%% global average

obs_d_m = nanmean(100*day_sa,2)';
obs_n_m = nanmean(100*night_sa,2)';

obs_d_std= nanstd(100*day_sa,1,2)';
obs_n_std = nanstd(100*night_sa,1,2)';

obs_d_q1 = quantile(100*day_sa,0.25,2)';
obs_d_q3 = quantile(100*day_sa,0.75,2)';
obs_n_q1 = quantile(100*night_sa,0.25,2)';
obs_n_q3 = quantile(100*night_sa,0.75,2)';

mod_d_m = nanmean(100*profile_total_d,2)';
mod_n_m = nanmean(100*profile_total_n,2)';

mod_d_std = nanstd(100*profile_total_d,1,2)';
mod_n_std = nanstd(100*profile_total_n,1,2)';

mod_d_q1 = quantile(100*profile_total_d,0.25,2)';
mod_d_q3 = quantile(100*profile_total_d,0.75,2)';
mod_n_q1 = quantile(100*profile_total_n,0.25,2)';
mod_n_q3 = quantile(100*profile_total_n,0.75,2)';

lab = {'obs','mod'};
x_max = 0;
for i=1:2
    eval(['profile_d_mean{i} =' lab{i} '_d_m;'])
    eval(['profile_n_mean{i} =' lab{i} '_n_m;'])

    
    eval(['profile_d_q1 =' lab{i} '_d_q1;'])
    eval(['profile_d_q3 =' lab{i} '_d_q3;'])
    eval(['profile_n_q1 =' lab{i} '_n_q1;'])
    eval(['profile_n_q3 =' lab{i} '_n_q3;'])
    
wmd_d_m(i) = nansum(profile_d_mean{i} .* depth/100);
wmd_n_m(i) = nansum(profile_n_mean{i} .* depth/100);


 x_fill_d{i} = [profile_d_q1,profile_d_q3(end:-1:1)];
 x_fill_n{i} = [profile_n_q1,profile_n_q3(end:-1:1)];


y_fill{i} = [depth, depth(end:-1:1)];

x_max = max(x_max,1.1 * max([x_fill_d{i} x_fill_n{i}]));

end

%% Figure 1-b

col1 = [255,51,153]/255;

fill(ax2, [-x_max, -x_max 0 0],[max(depth) min(depth) min(depth) max(depth)],'w')
hold(ax2,"on") 
fill(ax2,[0 0 x_max x_max],[max(depth) min(depth) min(depth) max(depth)],[0.5 0.5 0.5],'facealpha',0.2)

%obs
fill(ax2, -x_fill_d{1},y_fill{1},col1,'facealpha',0.3,'edgecolor','none')
fill(ax2, x_fill_n{1},y_fill{1},col1,'facealpha',0.3,'edgecolor','none')

hline(2) = plot(ax2, -profile_d_mean{1},depth,'color',col1,'linewidth',2,'DisplayName','Obs');
plot(ax2, profile_n_mean{1},depth,'color',col1,'linewidth',2)


%mod
fill(ax2, -x_fill_d{2},y_fill{2},'k','facealpha',0.2,'edgecolor','none')
fill(ax2, x_fill_n{2},y_fill{2},'k','facealpha',0.2,'edgecolor','none')

hline(1) = plot(ax2, -profile_d_mean{2},depth,'--k','linewidth',2,'DisplayName','Model');
plot(ax2, profile_n_mean{2},depth,'--k','linewidth',2)

xlim(ax2, [-7 7])
ylim(ax2, [-900 0])


x_tick = get(ax2, 'xtick');
y_tick = get(ax2, 'ytick');
x_tick_ = abs(x_tick);
x_tick_label = strsplit(num2str(x_tick_));
set(ax2, 'xtick',x_tick,'xticklabel',x_tick_label,'fontsize',14)
set(ax2, 'ytick',y_tick,'yticklabel',[])
xlabel(ax2, 'Vertical distribution (%)')

legend(hline,'box','off','location','northwest','fontsize',11)

% correlation
r = corrcoef(profile_d_mean{1},profile_d_mean{2});
r_square(1) = r(1,2)^2;
r = corrcoef(profile_n_mean{1},profile_n_mean{2});
r_square(2) = r(1,2)^2;

text(ax2,-4.5,-800,['R^2 = ' sprintf('%0.2f',r_square(1))],'fontsize',12,'color','k')
text(ax2,-4,-850,'Day','fontsize',12,'color','k')

text(ax2,2.25,-800,['R^2 = ' sprintf('%0.2f',r_square(2))],'fontsize',12,'color','k')
text(ax2,2.5,-850,'Night','fontsize',12,'color','k')


%% Figure 1-c

len = 0.04:0.02:0.22;
beta = -1;

d = importdata('contribution_com1.txt');
alpha1 = d.data;
alpha1 = alpha1(clstr);


d = importdata('contribution_com3.txt');
alpha3 = d.data;
alpha3 = alpha3(clstr);


alpha2 = 1 - alpha1 - alpha3;

load(['profile_by_com_'  prj ]) %profile_by_com

dn  = size(profile_by_com,1);
ng = size(profile_by_com,2);
nz = size(profile_by_com,3);
nx = size(profile_by_com,4);

for l=1:length(len)
      d = importdata(['sigma_' num2str(l) '.txt']);
      phi = d.data; phi = repmat(phi,1,dn,ng,nx); phi = permute(phi,[2 3 1 4]);
      toto(:,:,:,:,l) = phi * ((len(l)/len(1))^(4*beta));
end


profile_norm_ = squeeze(nansum(profile_by_com .* toto,5));
profile_norm_sum = squeeze(nansum(profile_norm_,3)); profile_norm_sum = repmat(profile_norm_sum,1,1,1,nz); profile_norm_sum = permute(profile_norm_sum,[1 2 4 3]);
profile_norm = profile_norm_./profile_norm_sum;


alpha1_ = repmat(alpha1',100,1); 
alpha3_ = repmat(alpha3',100,1); 
alpha2_ = 1 - alpha1_ - alpha3_;
alpha_ = cat(3,alpha1_,alpha2_,alpha3_); alpha_=permute(alpha_,[3 1 2]);


profile_com_d = 100*(alpha_.*squeeze(profile_norm(1,:,:,:)));
profile_com_n = 100*(alpha_.*squeeze(profile_norm(2,:,:,:)));

profile_com_d_m = squeeze(nanmean(profile_com_d,3));
profile_com_n_m = squeeze(nanmean(profile_com_n,3));

profile_com_d_std = squeeze(nanstd(profile_com_d,1,3));
profile_com_n_std = squeeze(nanstd(profile_com_n,1,3));

IC_d_min = profile_com_d_m - profile_com_d_std;
IC_d_max = profile_com_d_m + profile_com_d_std;
IC_n_min = profile_com_n_m - profile_com_n_std;
IC_n_max = profile_com_n_m + profile_com_n_std;

y_fill = [depth, depth(end:-1:1)];

x_max = 7;


col = {[255,215,0]/255 , [145,212,194]/255, [255,102,178]/255, [255,204,229]/255 , [175,238,238]/255 };

plot_with_shading= 0;


fill(ax3,[-x_max, -x_max 0 0],[max(depth) min(depth) min(depth) max(depth)],'w') 
hold(ax3,"on")
fill(ax3,[0 0 x_max x_max],[max(depth) min(depth) min(depth) max(depth)],[0.5 0.5 0.5],'facealpha',0.2)

leg = {'Epipelagic','Migrant','Resident'};


for i=1:size(IC_d_min,1)
    if(plot_with_shading)
        x_fill_d = [IC_d_min(i,:),IC_d_max(i,end:-1:1)];
        x_fill_n = [IC_n_min(i,:),IC_n_max(i,end:-1:1)];
        fill(ax3,x_fill_d,y_fill,col{i},'facealpha',0.2,'edgecolor','none')
        fill(ax3,-x_fill_n,y_fill,col{i},'facealpha',0.2,'edgecolor','none')
    end
    % day
    hLines(i) = plot(ax3,-profile_com_d_m(i,:),depth,'color',col{i},'linewidth',2,'DisplayName',leg{i}); %,'visible','off');
    y_axis = [depth,depth(end:-1:1)];

    x_axis = [profile_com_d_m(i,:),zeros(1,nz) ];
    a= fill(ax3, -x_axis, y_axis,col{i},'edgecolor',col{i},'linewidth',2,'facealpha',0.7);
        % night
    x_axis = [profile_com_n_m(i,:),zeros(1,nz) ];
    fill(ax3, x_axis, y_axis, col{i},'edgecolor',col{i},'linewidth',2,'facealpha',0.7)
    
    
end


hLines(4) = plot(ax3,-profile_d_mean{2},depth,'--','color','k','linewidth',1, 'DisplayName','Total');
plot(ax3,profile_n_mean{2},depth,'--','color','k','linewidth',1)


xlim(ax3,[-7 7])
ylim(ax3,[-900 0])

x_tick = get(ax3,'xtick');
x_tick_ = abs(x_tick);
x_tick_label = strsplit(num2str(x_tick_));
set(ax3,'xtick',x_tick,'xticklabel',x_tick_label,'fontsize',12)

set(ax3, 'ytick',-1000:200:0,'yticklabels', {'-1000    ', '-800    ','-600    ','-400    ', '-200    ','0    '}, 'fontsize',14)

xlabel(ax3, 'Vertical distribution (%)')

legend(hLines,'box','off','location','northwest','fontsize',11)

text(ax3,-4,-850,'Day','fontsize',12,'color','k')
text(ax3,2.5,-850,'Night','fontsize',12,'color','k')

title(ax1, "(a) Malaspina cruise track ",'FontSize',14,'FontWeight','normal')
title(ax2, "(b) All profiles: Obs. vs. Model ",'FontSize',14,'FontWeight','normal')
title(ax3, "(c) All profiles: Model communities ",'FontSize',14,'FontWeight','normal')


% 


