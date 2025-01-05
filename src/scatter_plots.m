
project = 'MALASPINA';


print_fig1 = 0;
print_fig2 = 0;

%scenario = [scenario suffix];


AB =  []; 

WMD_OBS_d = [];
WMD_MOD_d = [];
WMD_OBS_n = [];
WMD_MOD_n = [];
len = [0];
lat = [];
lon = [];
oxy_min = [];
    
%i_s = 4;

id_depth = 1:100;
    
pathway = 'data/';
   
load([pathway 'profile_total_' project '_' prj ])

load([pathway 'depth'])
id = find(depth<=1000);
depth = depth(id);
if(size(depth,1)>1)
    depth = depth';
end
depth = depth(id_depth);
% observations
d = importdata([pathway 'cluster_list_selected.txt']);
clstr = d.data;

load([pathway 'day_sa_m'])
load([pathway 'night_sa_m'])

day_sa = day_sa_m(id,clstr);
night_sa = night_sa_m(id,clstr);



ind_d = find(isnan(day_sa));
ind_n = find(isnan(night_sa));


load([pathway 'oxy_obs_m'])
oxy_ = nanmin(oxy_obs_m(1:100,clstr));


profile_total_d = squeeze(nanmean(profile_total(1,:,:,:),4));
profile_total_n = squeeze(nanmean(profile_total(2,:,:,:),4));

DEPTH = repmat(depth,2,1,length(clstr));
 
profile_total_d(ind_d) = nan;
profile_total_n(ind_n) = nan;

profile_total_d = profile_total_d(id_depth,:)./(nansum(profile_total_d(id_depth,:)));
profile_total_n = profile_total_n(id_depth,:)./(nansum(profile_total_n(id_depth,:)));


% model
WMD_mod_d = squeeze(nansum(profile_total_d.*squeeze(DEPTH(1,:,:))));
WMD_mod_n = squeeze(nansum(profile_total_n.*squeeze(DEPTH(1,:,:))));

WMD_mod = cat(3,WMD_mod_d,WMD_mod_n); WMD_mod = permute(WMD_mod,[3 1 2]);




oxy_min = [oxy_min,oxy_];

sa_d = day_sa(id_depth,:)./nansum(day_sa(id_depth,:));
sa_n = night_sa(id_depth,:)./nansum(night_sa(id_depth,:));
sa_obs = cat(3,sa_d,sa_n); sa_obs = permute(sa_obs,[3 1 2]);% 2*100*651


WMD_obs = squeeze(nansum(sa_obs.*squeeze(DEPTH(:,:,:,1)),2)); %2*651
% [~,ID] = max(sa_obs(:,21:100,:,:),[],2);
% D_max_obs = depth_d(squeeze(ID));

if(strcmp(project,'MALASPINA'))
    WMD_obs(:,AB) = [];
end

WMD_OBS_d = [WMD_OBS_d,WMD_obs(1,:)];
WMD_MOD_d = [WMD_MOD_d,WMD_mod(1,:)];
WMD_OBS_n = [WMD_OBS_n,WMD_obs(2,:)];
WMD_MOD_n = [WMD_MOD_n,WMD_mod(2,:)];
len = [len, len + size(WMD_obs,2)];


clearvars sa_obs ind_d ind_n %-except PRJ prj scenario p_value


%%
wmd_obs = [];
wmd_mod = [];

thr_value = -50;

load('data/id_')
load('data/id_over_est')
load('data/id_under_est')

col_id = {[75 0 130]/225,[0 1 0],[1 0 0]};

%%



%ax1 = subplot_tight(1,2,1,[0.2, 0.05]);
if(plot_fig3_d)
    fig = figure();
   ax1 = axes(fig);
else
    fig = figure('Renderer', 'painters', 'Position', [10 30 1000 400]);
    ax1 = subplot(1,2,1);
    ax2 = subplot(1,2,2);
    set(ax2,'box' , 'on');
end

set(ax1,'box' , 'on');

%ax2 = subplot_tight(1,2,2,[0.2, 0.05]);

%%

if(plot_fig2_a_b || plot_fig3_d || plot_fig4_a_b)

    

    colorOrder = get(gca, 'ColorOrder');
    colorOrder(4,:) = [1,0,0];
    marker_size= 40;

    plot(ax1,[-900 0],[-900 0],'k','linewidth',2)
    hold(ax1,'on')
    plot(ax1,[-900-thr_value 0],[-900 thr_value],'--k','linewidth',1)
    plot(ax1,[-900 thr_value],[-900-thr_value 0],'--k','linewidth',1)

        wmd_obs_ = WMD_OBS_d(len(1)+1:len(2));
        wmd_mod_ = WMD_MOD_d(len(1)+1:len(2));
        oxy_min_ = oxy_min(len(1)+1:len(2));
    
    
        wmd_diff = wmd_obs_ - wmd_mod_;
    
        if(strcmp(prj , 'three_groups_without_oxy'))
        
            id_over_est = find(wmd_diff<=thr_value);
            id_under_est = find(wmd_diff>=-thr_value);
        
             id=1:length(wmd_diff);
        
            id_ = id(~ismember(id,id_over_est) & ~ismember(id,id_under_est));

        end
              
  
        scatter(ax1,-wmd_obs_(id_),-wmd_mod_(id_),marker_size,[75 0 130]/225,'o','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);%,'facealpha',0.5);
  
        scatter(ax1,-wmd_obs_(id_over_est),-wmd_mod_(id_over_est),marker_size,[0 86 27]/225,'o','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.8); 

        scatter(ax1,-wmd_obs_(id_under_est),-wmd_mod_(id_under_est),marker_size,colorOrder(4,:),'o','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);

%  
        wmd_obs = [wmd_obs,wmd_obs_(id)];
        wmd_mod = [wmd_mod,wmd_mod_(id)];

        xlim(ax1,[-900 0]);
        ylim(ax1,[-900 0])

        xlabel(ax1,'Observed WMD (m)')
        


        [R,p_v] = corrcoef(wmd_obs,wmd_mod);
        R2 = (round(R(1,2)*100)/100)^2;
        R1 = (round(R(1,2)*100)/100);


        if(plot_fig3_d)
                text(ax1,-600,-800,['R^2 = 0.71 -> ' sprintf('%0.2f',R2)],'fontsize',12)
                text(ax1,-500,-850,"Day",'fontsize',12)
                title(ax1,"(d) All WMD: O2 impact",'FontWeight','normal','fontsize',12)
        elseif(plot_fig4_a_b)
                text(ax1,-400,-800,['R^2 = 0.71 -> ' sprintf('%0.2f',R2)],'fontsize',12)
                title(ax1,"(a) Daytime WMDs")
        elseif(plot_fig2_a_b)

                text(ax1,-300,-800,['R^2 = ' sprintf('%0.2f',R2)],'fontsize',13)
                title(ax1,"(a) Daytime WMDs",'fontsize',12,'FontWeight','normal')
        end
        set(ax1,'YTick',-900:100:0,'yTickLabel',{'','','','','','','','','',''});
        set(ax1,'XTick',-900:100:0,'xTickLabel',{'','-800','','-600','','-400','','-200','','0'},'XTickLabelRotation',0,'fontsize',12);

        ylabel(ax1, 'Predicted WMD (m)', 'fontsize',14)

        % for legend:
        scatter(ax1,-850,-90,150,col_id{1},'o','filled','markeredgecolor','k')
        text(ax1, -800, -90, 'GP stations','color',col_id{1},'fontsize',12)
        scatter(ax1,-850,-170,150,col_id{2},'o','filled','markeredgecolor','k')
        text(ax1, -800, -170, 'SP stations','color',col_id{2},'fontsize',12)
        scatter(ax1,-850,-250,150,col_id{3},'o','filled','markeredgecolor','k')
        text(ax1, -800, -250, 'DP stations','color',col_id{3},'fontsize',12)
end


%%

if((plot_fig2_a_b || plot_fig4_a_b) && ~plot_fig3_d )

    wmd_obs = [];
    wmd_mod = [];

    wmd_obs_ = WMD_OBS_n(len(1)+1:len(2));
    wmd_mod_ = WMD_MOD_n(len(1)+1:len(2));
    oxy_min_ = oxy_min(len(1)+1:len(2));

  
    

    plot(ax2,[-900 0],[-900 0],'k','linewidth',2)
    hold(ax2,'on')
    plot(ax2,[-900-thr_value 0],[-900 thr_value],'--k','linewidth',1)
    plot(ax2,[-900 thr_value],[-900-thr_value 0],'--k','linewidth',1)
  

    scatter(ax2,-wmd_obs_(id_),-wmd_mod_(id_),marker_size,[75 0 130]/225,'o','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);%,'facealpha',0.5);
    scatter(ax2,-wmd_obs_(id_over_est),-wmd_mod_(id_over_est),marker_size,[0 86 27]/225,'o','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);  

    scatter(ax2,-wmd_obs_(id_under_est),-wmd_mod_(id_under_est),marker_size,colorOrder(4,:),'o','filled','MarkerEdgeColor','k','MarkerFaceAlpha',0.8);

  
    wmd_obs = [wmd_obs,wmd_obs_(id)];
    wmd_mod = [wmd_mod,wmd_mod_(id)];

    xlim(ax2,[-900 0])
    ylim(ax2,[-900 0])

    set(ax2,'color',[240 240 240]/255);
    xlabel(ax2, 'Observed WMD (m)')
    %ylabel(ax2, 'Predicted WMD (m)')
    set(ax2,'YTick',-900:100:0,'yTickLabel',{'','-800  ','' ,'-600  ','','-400  ','','-200  ','','0  '}, 'fontsize',12);
    set(ax2,'XTick',-900:100:0,'xTickLabel',{'','-800  ','' ,'-600  ','','-400  ','','-200  ','','0  '}, 'XTickLabelRotation',0,'fontsize',12)
    
    [R,p_v] = corrcoef(wmd_obs,wmd_mod);
    R2 = (round(R(1,2)*100)/100)^2;
    R1 = (round(R(1,2)*100)/100);
    

    if(plot_fig2_a_b)
          text(-300,-800,['R^2 = ' sprintf('%0.2f',R2)],'fontsize',13)
          title("(b) Nighttime WMDs",'fontsize',12,'FontWeight','normal')
    elseif(plot_fig4_a_b)
          text(-400,-800,['R^2 = 0.94 -> ' sprintf('%0.2f',R2)],'fontsize',12)
          title("(b) Nighttime WMDs",'fontsize',12,'FontWeight','normal')
     end


end


    