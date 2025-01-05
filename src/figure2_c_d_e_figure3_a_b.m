%clearvars;

%% PARAMETERS
project = 'MALASPINA';
%prj = 'three_groups_without_oxy';
print_fig2_a = 0;
print_fig2_b = 0;
print_fig2_c = 0;

print_fig3_a = 0;
print_fig3_b = 0;

print_fig2 = {'print_fig2_a','print_fig2_b','print_fig2_c'};
pathway = 'data/';

load([pathway 'id_'])
load([pathway 'id_over_est'])
load([pathway 'id_under_est'])

load([pathway 'depth'])
id = find(depth<=1000);
if(size(depth,1)>1)
    depth = depth';
end
id_depth = 1:100;
depth = depth(id_depth);

d = importdata([pathway 'cluster_list_selected.txt']);
clstr = d.data;

load([pathway 'day_sa_m'])
load([pathway 'night_sa_m'])


day_sa = day_sa_m(id,clstr);
night_sa = night_sa_m(id,clstr);



load([pathway 'par_obs_m'])

par_d = par_obs_m(:,clstr);

load([pathway 'depth'])
id_dep = find(depth<=1000);
depth = depth(id_dep);
depth_ = -depth;


%% 

profile_d_mean_reg{1} = nanmean(day_sa(:,id_),2);
profile_d_mean_reg{2} = nanmean(day_sa(:,id_over_est),2);
profile_d_mean_reg{3} = nanmean(day_sa(:,id_under_est),2);
profile_d_std_reg{1} = nanstd(day_sa(:,id_),1,2);
profile_d_std_reg{2} = nanstd(day_sa(:,id_over_est),1,2);
profile_d_std_reg{3} = nanstd(day_sa(:,id_under_est),1,2);

profile_d_q1_reg{1} = quantile(day_sa(:,id_),0.25,2);
profile_d_q1_reg{2} = quantile(day_sa(:,id_over_est),0.25,2);
profile_d_q1_reg{3} = quantile(day_sa(:,id_under_est),0.25,2);
profile_d_q3_reg{1} = quantile(day_sa(:,id_),0.75,2);
profile_d_q3_reg{2} = quantile(day_sa(:,id_over_est),0.75,2);
profile_d_q3_reg{3} = quantile(day_sa(:,id_under_est),0.75,2);



N(1) = length(id_);
N(2) = length(id_over_est);
N(3) = length(id_under_est);


profile_n_mean_reg{1} = nanmean(night_sa(:,id_),2);
profile_n_mean_reg{2} = nanmean(night_sa(:,id_over_est),2);
profile_n_mean_reg{3} = nanmean(night_sa(:,id_under_est),2);

profile_n_std_reg{1} = nanstd(night_sa(:,id_),1,2);
profile_n_std_reg{2} = nanstd(night_sa(:,id_over_est),1,2);
profile_n_std_reg{3} = nanstd(night_sa(:,id_under_est),1,2);

profile_n_q1_reg{1} = quantile(night_sa(:,id_),0.25,2);
profile_n_q1_reg{2} = quantile(night_sa(:,id_over_est),0.25,2);
profile_n_q1_reg{3} = quantile(night_sa(:,id_under_est),0.25,2);
profile_n_q3_reg{1} = quantile(night_sa(:,id_),0.75,2);
profile_n_q3_reg{2} = quantile(night_sa(:,id_over_est),0.75,2);
profile_n_q3_reg{3} = quantile(night_sa(:,id_under_est),0.75,2);


par_d_mean_reg{1} = log10(nanmean(par_d(:,id_),2));
par_d_mean_reg{2} = log10(nanmean(0.5*par_d(:,id_over_est),2));
par_d_mean_reg{3} = log10(nanmean(par_d(:,id_under_est),2));



col_id = {[75 0 130]/225,[0 86 27]/225,[1 0 0]};
%col_id = {[75 0 130]/225,[0 40 25]/225,[1 0 0]};
%%
if(plot_fig3_c)
    fig2 = figure('visible','off');
    ax_fig2 = axes();
end

if(plot_fig3_a_b)
    
    
    fig = figure('Renderer', 'painters', 'Position', [10 30 1200 500]);
    
    

    ax1 = subplot(2,4,[1,2,5,6]);
    ax2 = subplot(2,4,[3,4,7,8]);

    leg_lab = {'GP','SP','DP'};

%if(plot_fig3_a)

    %figure;

    for i=1:3
        par_v = par_d_mean_reg{i};
        IC_d_min = profile_d_q1_reg{i};
        IC_d_max = profile_d_q3_reg{i};
        x_fill = [IC_d_min',IC_d_max(end:-1:1)'];
        yfill = [par_v',par_v(end:-1:1)'];

    
        f(i) = fill(ax1, 100*x_fill,yfill,col_id{i},'facealpha',0.3,'edgecolor','none','DisplayName',leg_lab{i});
        hold(ax1,'on')
        plot(ax1, 100*profile_d_mean_reg{i},par_v,'color',col_id{i},'linewidth',2)
    
        % les stats
        V = profile_d_mean_reg{i}';
        [MODE(i),Qrt{i}] = estimate_statistic_metrics(V(21:100),par_v(21:100));

        
    
    end
    xlim(ax1,[0 13])
    ylim(ax1,[-10 2])
    ylabel(ax1,'log_{10}(PAR)')
    xlabel(ax1,'Vertical distribution (%)')

    title(ax1,'(a) Acoustic profiles vs. PAR','fontsize',14,'fontweight','normal')
    
    set(ax1,'ytick',ax1.get('ytick'),'yticklabel',[],'fontsize',12)

    lgd = legend(ax1,f,'fontsize',16,'location','northeast','box','off');

end
%% env variables

load([pathway 'oxy_obs_m'])
oxy_obs_m_ = oxy_obs_m(:,clstr);

oxy_mean_reg{1} = nanmean(oxy_obs_m_(:,id_),2);
oxy_mean_reg{2} = nanmean(oxy_obs_m_(:,id_over_est),2);
oxy_mean_reg{3} = nanmean(oxy_obs_m_(:,id_under_est),2);
oxy_std_reg{1} = nanstd(oxy_obs_m_(:,id_),1,2);
oxy_std_reg{2} = nanstd(oxy_obs_m_(:,id_over_est),1,2);
oxy_std_reg{3} = nanstd(oxy_obs_m_(:,id_under_est),1,2);

oxy_q1_reg{1} = quantile(oxy_obs_m_(:,id_),0.25,2);
oxy_q1_reg{2} = quantile(oxy_obs_m_(:,id_over_est),0.25,2);
oxy_q1_reg{3} = quantile(oxy_obs_m_(:,id_under_est),0.25,2);
oxy_q3_reg{1} = quantile(oxy_obs_m_(:,id_),0.75,2);
oxy_q3_reg{2} = quantile(oxy_obs_m_(:,id_over_est),0.75,2);
oxy_q3_reg{3} = quantile(oxy_obs_m_(:,id_under_est),0.75,2);
oxy_median_reg{1} = quantile(oxy_obs_m_(:,id_),0.5,2);
oxy_median_reg{2} = quantile(oxy_obs_m_(:,id_over_est),0.5,2);
oxy_median_reg{3} = quantile(oxy_obs_m_(:,id_under_est),0.5,2);


if(plot_fig3_a_b)

    %subplot(121)
    for i=1:3
        par_v = par_d_mean_reg{i};
        IC_d_min = oxy_q1_reg{i};
        IC_d_max = oxy_q3_reg{i};
        x_fill = [IC_d_min',IC_d_max(end:-1:1)'];
        yfill = [par_v',par_v(end:-1:1)'];
        
    
        fill(ax2, x_fill,yfill,col_id{i},'facealpha',0.3,'edgecolor','none')
        hold(ax2,'on')
        plot(ax2, oxy_mean_reg{i},par_v,'color',col_id{i},'linewidth',2)   
    
    end
    xlim(ax2,[0 6])
    ylim(ax2,[-10 2])
    xlabel(ax2,'Oxygen (ml l^{-1})')
    set(ax2,'fontsize',12)
    title(ax2,'(b) Oxygen vs. PAR','fontsize',14,'fontweight','normal')
    

   
    

end

%%



%%

id_over_est = [84:95];

%%

load([pathway 'profile_total_' project '_' prj ])
profile_total_d = squeeze(nanmean(profile_total(1,:,:,:),4));
profile_total_n = squeeze(nanmean(profile_total(2,:,:,:),4));

if(add_simple_model)
    load([pathway 'profile_total_' project '_three_groups_without_oxy'])
    profile_total_simple_d = squeeze(nanmean(profile_total(1,:,:,:),4));
    profile_total_simple_n = squeeze(nanmean(profile_total(2,:,:,:),4));

end


profile_d_mean_reg{1} = nanmean(day_sa(:,id_),2);
profile_d_mean_reg{2} = nanmean(day_sa(:,id_over_est),2);
profile_d_mean_reg{3} = nanmean(day_sa(:,id_under_est),2);
profile_d_std_reg{1} = nanstd(day_sa(:,id_),1,2);
profile_d_std_reg{2} = nanstd(day_sa(:,id_over_est),1,2);
profile_d_std_reg{3} = nanstd(day_sa(:,id_under_est),1,2);

profile_d_q1_reg{1} = quantile(day_sa(:,id_),0.25,2);
profile_d_q1_reg{2} = quantile(day_sa(:,id_over_est),0.25,2);
profile_d_q1_reg{3} = quantile(day_sa(:,id_under_est),0.25,2);
profile_d_q3_reg{1} = quantile(day_sa(:,id_),0.75,2);
profile_d_q3_reg{2} = quantile(day_sa(:,id_over_est),0.75,2);
profile_d_q3_reg{3} = quantile(day_sa(:,id_under_est),0.75,2);


profile_n_mean_reg{1} = nanmean(night_sa(:,id_),2);
profile_n_mean_reg{2} = nanmean(night_sa(:,id_over_est),2);
profile_n_mean_reg{3} = nanmean(night_sa(:,id_under_est),2);

profile_n_std_reg{1} = nanstd(night_sa(:,id_),1,2);
profile_n_std_reg{2} = nanstd(night_sa(:,id_over_est),1,2);
profile_n_std_reg{3} = nanstd(night_sa(:,id_under_est),1,2);

profile_n_q1_reg{1} = quantile(night_sa(:,id_),0.25,2);
profile_n_q1_reg{2} = quantile(night_sa(:,id_over_est),0.25,2);
profile_n_q1_reg{3} = quantile(night_sa(:,id_under_est),0.25,2);
profile_n_q3_reg{1} = quantile(night_sa(:,id_),0.75,2);
profile_n_q3_reg{2} = quantile(night_sa(:,id_over_est),0.75,2);
profile_n_q3_reg{3} = quantile(night_sa(:,id_under_est),0.75,2);



  label = {'','over_est','under_est'};
  label_ = {'profile','mod'};
  fig_suffix = {'a','b','c'};



    r_coef = [0.89, 0.56, 0.51 ; 0.87, 0.74, 0.77];

    x_max_ = [7, 20, 6];

    yfill = [depth_,depth_(end:-1:1)];

    title_lab = {'(c) GP profiles','(d) SP profiles', '(e) DP profiles'};

if(plot_fig2_c_d_e)

    fig = figure('Renderer', 'painters', 'Position', [10 30 1200 500]);

end

for i_reg = 1:length(label)

    if(plot_fig2_c_d_e)
        eval(['ax = subplot_tight(1,3,' num2str(i_reg) ',[0.2, 0.035]);']) %[0.08, 0.05]

        set(ax, 'box' , 'on');

    else
        ax = gca;
    end

    eval(['id_reg = id_' label{i_reg} ';'])

    if(add_simple_model)
        mod_simple_d_mean_reg{i_reg} = nanmean(profile_total_simple_d(:,id_reg),2);
        mod_simple_n_mean_reg{i_reg} = nanmean(profile_total_simple_n(:,id_reg),2);
    end
    
    
    mod_d_mean_reg{i_reg} = nanmean(profile_total_d(:,id_reg),2);
    mod_d_std_reg{i_reg} = nanstd(profile_total_d(:,id_reg),1,2);
    mod_n_mean_reg{i_reg} = nanmean(profile_total_n(:,id_reg),2);
    mod_n_std_reg{i_reg} = nanstd(profile_total_n(:,id_reg),1,2);
    
    mod_d_q1_reg{i_reg} = quantile(profile_total_d(:,id_reg),0.25,2);
    mod_d_q3_reg{i_reg} = quantile(profile_total_d(:,id_reg),0.75,2);
    mod_n_q1_reg{i_reg} = quantile(profile_total_n(:,id_reg),0.25,2);
    mod_n_q3_reg{i_reg} = quantile(profile_total_n(:,id_reg),0.75,2);
    
    for j=1:length(label_)
        eval(['IC_d_min =' label_{j} '_d_q1_reg{i_reg};']) %' label_{j} '_d_mean_reg{i_reg} - label_{j} '_d_std_reg{i_reg}/sqrt(N(i_reg));']) %
        eval(['IC_d_max =' label_{j} '_d_q3_reg{i_reg};']) %' label_{j} '_d_mean_reg{i_reg} + label_{j} '_d_std_reg{i_reg}/sqrt(N(i_reg));']) %
        eval(['IC_n_min =' label_{j}  '_n_q1_reg{i_reg};']) %' label_{j} '_n_mean_reg{i_reg} -label_{j} '_n_std_reg{i_reg}/sqrt(N(i_reg));']) %
        eval(['IC_n_max =' label_{j} '_n_q3_reg{i_reg};']) %' label_{j} '_n_mean_reg{i_reg} + label_{j} '_n_std_reg{i_reg}/sqrt(N(i_reg));']) %
        
        x_fill_d = [IC_d_min',IC_d_max(end:-1:1)'];
        x_fill_n = [IC_n_min',IC_n_max(end:-1:1)'];
        
        eval(['x_fill_d_' label_{j}  '= x_fill_d;'])
        eval(['x_fill_n_' label_{j}  '= x_fill_n;'])
        
        %x_max_(j) = max([x_fill_d,x_fill_n]);
    end
        %x_max= 1.2*100*max(x_max_)

        x_max = x_max_(i_reg);%


        % correlation coefficient between predicted and observed profiles    
        r = corrcoef(mod_d_mean_reg{i_reg},profile_d_mean_reg{i_reg});
        corr_coef_d(i_reg) = r(1,2)^2;
   
    
        r = corrcoef(mod_n_mean_reg{i_reg},profile_n_mean_reg{i_reg});
        corr_coef_n(i_reg) = r(1,2)^2;
     
  if(plot_fig2_c_d_e)

     
        
        fill(ax,[0 0 x_max x_max],[max(depth_) min(depth_) min(depth_) max(depth_)],[0.5 0.5 0.5],'facealpha',0.2)
        hold on
        fill(ax,[-x_max, -x_max 0 0],[max(depth_) min(depth_) min(depth_) max(depth_)],'w') %,'facealpha',0.05)
     
        hl(1) = fill(ax, -100*x_fill_d_mod,yfill,'k','facealpha',0.3,'edgecolor','none','DisplayName','Model');
        hold(ax,'on')
        fill(ax, 100*x_fill_n_mod,yfill,'k','facealpha',0.3,'edgecolor','none')

        hl(2) =fill(ax, -100*x_fill_d_profile,yfill,col_id{i_reg},'facealpha',0.3,'edgecolor','none','DisplayName','Obs');
        fill(ax, 100*x_fill_n_profile,yfill,col_id{i_reg},'facealpha',0.3,'edgecolor','none')
        
        
        plot(ax, -100*mod_d_mean_reg{i_reg},depth_,'--k','linewidth',1.5)
        plot(ax, 100*mod_n_mean_reg{i_reg},depth_,'--k','linewidth',1.5)
        plot(ax, -100*profile_d_mean_reg{i_reg},depth_,'color',col_id{i_reg},'linewidth',2)
        plot(ax, 100*profile_n_mean_reg{i_reg},depth_,'color',col_id{i_reg},'linewidth',2)

        
        if(add_simple_model)
           hl(3) = plot(ax, -100*mod_simple_d_mean_reg{i_reg},depth_,'--k','linewidth',1,'DisplayName','Simple model');
           plot(ax, 100*mod_simple_n_mean_reg{i_reg},depth_,'--k','linewidth',1) 
        end

        [hl,lgd] = legend(hl,'location','NorthWest','box','off','fontsize',12);
        if(plot_fig4_c_d_e)
            lgd(4).FaceAlpha = 0.3;
            lgd(5).FaceAlpha = 0.3;
        else
            lgd(3).FaceAlpha = 0.3;
            lgd(4).FaceAlpha = 0.3;
        end
    
        xlim(ax, [-x_max x_max])%
        ylim(ax, [-900 0])
        x_tick = get(ax,'xtick');
        x_tick_ = abs(x_tick);
        x_tick_label = strsplit(num2str(x_tick_));
        set(ax,'xtick',x_tick,'xticklabel',x_tick_label) %,'fontweight','bold','fontsize',12)
        title(ax,title_lab{i_reg}, 'FontSize',14,'FontWeight','normal')
    
        
        xlabel(ax, 'Vertical distribution (%)',"FontSize",14)
        %grid minor
        %set(ax,'fontsize',9,'fontweight','bold')
        %set(gcf, 'InvertHardCopy', 'off');

        if(i_reg==1)
            set(ax,'ytick',-900:100:0,'yticklabel',[])
            ylabel(ax, 'Depth (m)',"FontSize",14)
        else
            set(ax,'ytick',-900:100:0,'yticklabel',{'','-800','','-600','','-400','','-200','','0'},'fontsize',11)
        
        end
   
        

        if(plot_fig4_c_d_e)
            text(ax, 0.8*x_max/2,-800,'Night','fontsize',12)
            text(ax, 0.5*x_max/2,-850,['R^2 = ' num2str(r_coef(2,i_reg)) ' -> ' sprintf('%0.2f',corr_coef_n(i_reg))],'fontsize',10)

            text(ax, -1.2*x_max/2,-800,'Day','fontsize',12)
            text(ax, -1.75*x_max/2,-850,['R^2 = ' num2str(r_coef(1,i_reg)) ' -> ' sprintf('%0.2f',corr_coef_d(i_reg))],'fontsize',10)
        else    

            text(ax, 0.6*x_max/2,-800,['R^2 = ' sprintf('%0.2f',corr_coef_n(i_reg))],'fontsize',12)
            text(ax, 0.8*x_max/2,-850,'Night','fontsize',12)

            text(ax, -1.4*x_max/2,-800,['R^2 = ' sprintf('%0.2f',corr_coef_d(i_reg))],'fontsize',12)
            text(ax, -1.2*x_max/2,-850,'Day','fontsize',12)
        end
 

        eval(['print_fig =' print_fig2{i_reg} ';']);

        %savefig(['figures/figure2_' fig_suffix{i_reg}])

        

        if(print_fig)
            set(gcf,'renderer','Painters')
            fig_name = print_fig2{i_reg};
            exportgraphics(gcf,[fig_name '.eps'],'BackgroundColor','none','ContentType','vector')
            print -depsc -tiff -r300 -painters fig_name ;


            print(fig_name,'-dpng')

        end

  end

  if(plot_fig3_c && i_reg==2)

        fill(ax_fig2,[0 0 x_max x_max],[max(depth_) min(depth_) min(depth_) max(depth_)],[0.5 0.5 0.5],'facealpha',0.2)
        hold on
        fill(ax_fig2,[-x_max, -x_max 0 0],[max(depth_) min(depth_) min(depth_) max(depth_)],'w') %,'facealpha',0.05)
     
        hl(1) = fill(ax_fig2, -100*x_fill_d_mod,yfill,'k','facealpha',0.3,'edgecolor','none','DisplayName','Model');
        hold(ax_fig2,'on')
        fill(ax_fig2, 100*x_fill_n_mod,yfill,'k','facealpha',0.3,'edgecolor','none')

        hl(2) =fill(ax_fig2, -100*x_fill_d_profile,yfill,col_id{i_reg},'facealpha',0.3,'edgecolor','none','DisplayName','Obs');
        fill(ax_fig2, 100*x_fill_n_profile,yfill,col_id{i_reg},'facealpha',0.3,'edgecolor','none')
        
        
        plot(ax_fig2, -100*mod_d_mean_reg{i_reg},depth_,'--k','linewidth',0.5)
        plot(ax_fig2, 100*mod_n_mean_reg{i_reg},depth_,'--k','linewidth',0.5)
        plot(ax_fig2, -100*profile_d_mean_reg{i_reg},depth_,'color',col_id{i_reg},'linewidth',2)
        plot(ax_fig2, 100*profile_n_mean_reg{i_reg},depth_,'color',col_id{i_reg},'linewidth',2) 

        text(ax_fig2, 0.8*x_max/2,-850,'Night','fontsize',12)
        text(ax_fig2, 0.3*x_max/2,-800,['R^2 = 0.74 -> ' sprintf('%0.2f',corr_coef_n(i_reg))],'fontsize',12)

        text(ax_fig2, -1.2*x_max/2,-850,'Day','fontsize',12)
        text(ax_fig2, -1.75*x_max/2,-800,['R^2 = 0.56 -> ' sprintf('%0.2f',corr_coef_d(i_reg))],'fontsize',12)

        if(add_simple_model)
           plot(ax, -100*mod_simple_d_mean_reg{i_reg},depth_,'k','linewidth',1);
           plot(ax, 100*mod_simple_n_mean_reg{i_reg},depth_,'--k','linewidth',1) 
        end



        fig2.Visible = 'on';

        xlim(ax_fig2, [-x_max x_max])%
        ylim(ax_fig2, [-900 0])
        x_tick = get(ax_fig2,'xtick');
        x_tick_ = abs(x_tick);
        x_tick_label = strsplit(num2str(x_tick_));
        set(ax_fig2,'xtick',x_tick,'xticklabel',x_tick_label,'fontsize',12) %,'fontweight','bold',)
    
        set(ax_fig2,'ytick',-900:100:0,'yticklabel',{}) %,{'','-800','','-600','','-400','','-200','','0'},'fontsize',12)

        ylabel(ax_fig2,'Depth (m)','FontSize',14)
        
        xlabel(ax_fig2, 'Vertical distribution (%)',"FontSize",14)
        
        title(ax_fig2,'(c) SP: O2 impact','fontsize',14,'fontweight','normal')

        [hl,lgd] = legend(hl,'location','NorthWest','box','off','fontsize',14);
        lgd(3).FaceAlpha = 0.3;
        lgd(4).FaceAlpha = 0.3;

  end



end



 %%
 ox = oxy_mean_reg{2};
 
 V = mod_d_mean_reg{2}; % le profile prédit moyen dans la région OMZ
 V_cum = cumsum(V);
 [val_25, id_val_25] = min(abs(V_cum - 0.25));
 [val_50, id_val_50] = min(abs(V_cum - 0.50));
 [val_75, id_val_75] = min(abs(V_cum - 0.75));
 
 dep_25 = depth_(id_val_25);
 dep_50 = depth_(id_val_50);
 dep_75 = depth_(id_val_75);
 
 [val_max, id_val_max] = max(V);
 dep_max = depth_(id_val_max);
 oxy_depth__max = ox(id_val_max);
 
 v_obs = profile_d_mean_reg{2};
 V_cum_obs = cumsum(v_obs);
 [val_o_25, id_val_o_25] = min(abs(V_cum_obs - 0.25));
 [val_o_50, id_val_o_50] = min(abs(V_cum_obs - 0.50));
 [val_o_75, id_val_o_75] = min(abs(V_cum_obs - 0.75));
 
 dep_o_25 = depth_(id_val_o_25);
 dep_o_50 = depth_(id_val_o_50);
 dep_o_75 = depth_(id_val_o_75);
 
 [val_obs_max, id_val_obs_max] = max(v_obs);
 dep_max_obs = depth_(id_val_obs_max);
 oxy_depth__max_obs = ox(id_val_obs_max);
 
 [r_,p_] = corrcoef(V,v_obs); r = r_(1,2); p_val = p_(1,2);
