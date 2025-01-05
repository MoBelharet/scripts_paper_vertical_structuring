%project = 'MALASPINA';
%PRJ = {'three_groups_without_oxy','three_groups_without_oxy_'};%;
%prj = 'five_groups_bis';%'three_groups_without_oxy';
%pathway = ['/media/belharet/HD_belharet/Optimization_admb/outputs/' project '/transect/'];
%pathway_ = ['/media/belharet/HD_belharet/Optimization_admb/data/' project '/'];  

%d = importdata([pathway_ 'cluster_list_selected.txt']);
%clstr = d.data;

%id_depth = 1:100;

%load([pathway 'depth'])
%depth = depth(id_depth);

% observations
%load([pathway 'day_sa_m'])
%load([pathway 'night_sa_m'])

%day_sa = day_sa_m(id_depth,clstr);
%night_sa = night_sa_m(id_depth,clstr);


PRJ = {'five_groups_with_oxy','three_groups_without_oxy'};
suffix ={'_','_',''}; %_no_oxy

load('id_')
load('id_over_est')
id_over_est = [84:95];
load('id_under_est')


for i_prj = 1:length(PRJ)
    prj = PRJ{i_prj};
    load(['profile_total_MALASPINA_' prj ])
    
    profile_total_d = squeeze(profile_total(1,:,:));
    profile_total_n = squeeze(profile_total(2,:,:));
    
    %profile_total_d = profile_total_d_(id_depth,:)./(nansum(profile_total_d_(id_depth,:)));
    %profile_total_n = profile_total_n_(id_depth,:)./(nansum(profile_total_n_(id_depth,:)));
    
    for i_st = 1:length(clstr)
        
        mod_d = profile_total_d(:,i_st);
        mod_n = profile_total_n(:,i_st);
        obs_d = day_sa(:,i_st);
        obs_n = night_sa(:,i_st);
        
        r_d = corrcoef(mod_d,obs_d);
        r_square_d(i_prj,i_st) = r_d(1,2)^2 ; 
        
        r_n = corrcoef(mod_n,obs_n);
        r_square_n(i_prj,i_st) = r_n(1,2)^2 ;
        
        rms_d(i_prj,i_st) = 100*sqrt(nanmean(  ( ( mod_d-mean(mod_d) )-(obs_d-mean(obs_d)  )).^2));
        rms_n(i_prj,i_st) = 100*sqrt(nanmean(  ( ( mod_n-mean(mod_n) )-(obs_n-mean(obs_n)  )).^2));

        p = polyfit(mod_d,obs_d,1);
        slope_d(i_prj,i_st) = p(1);

        p = polyfit(mod_n,obs_n,1);
        slope_n(i_prj,i_st) = p(1);
        
    end
    
    
    
end