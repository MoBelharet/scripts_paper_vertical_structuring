%clearvars;

addpath('data/')
%%

titre = {'GP group', 'SP Group','DP group'};
couleur = {[75 0 130]/225,[0 86 27]/225,[1 0 0]};

legend_lab = {'Deep Res','Deep Mig','Shallow Res','Shallow Mig','Epi','Location','north'};

%%

fig = figure('Renderer', 'painters', 'Position', [20 30 1200 600]);


for i=1:3

    eval(['ax' num2str(i) '= subplot_tight(1,3,' num2str(i) ',[0.2, 0.05]);']) %[0.08, 0.05]

    eval(['set(ax' num2str(i) ',' '''box''' ',' '''on''' ');'])

    eval(['ax = ax' num2str(i) ';'])


end

%%
col1 = [255,215,0]/255;
col2 = [145,212,194]/255;
col3 = [255,102,178]/255;
col5 = [255,204,229]/255;
col4 = [175,238,238]/255;

%%
prj = 'five_groups_with_oxy';%'three_groups_without_oxy';%

load('id_')
id_over_est = 84:95;
load('id_under_est')

d = importdata('cluster_list_selected.txt');
clstr_list = d.data;

[clstr, id_sort] = sort(clstr_list,'ascend');

load('day_sa_m_NoNorm')
load('night_sa_m_NoNorm')
sa_d_t = nansum(day_sa_m(:,clstr)); % integrated over z
sa_n_t = nansum(night_sa_m(:,clstr));


dn = 2;
ng = 5;
nx = length(clstr);
nz = 100;
len = 0.04:0.02:0.22;
beta = -1;

acoustic_day = repmat(sa_d_t,ng,1,nz); acoustic_day = permute(acoustic_day,[1,3,2]);
acoustic_night = repmat(sa_n_t,ng,1,nz); acoustic_night = permute(acoustic_night,[1,3,2]);

d = importdata('contribution_com1.txt');
alpha1 = d.data;
alpha1 = alpha1(clstr);



d = importdata('contribution_com3.txt');
alpha3 = d.data;
alpha3 = alpha3(clstr);


alpha2 = 1 - alpha1 - alpha3;


if(ng==5)
    load(['P_com4_' prj]) % p
    load(['P_com5_' prj]) % q

    if(strcmp(prj,'five_groups_bis') )
            q(85) = 0.9998;
            q(91) = 0.3;
            q(92) = 0.5;
            q(93)=  0.9;
            q(95) = 0.6;
            q(98) = 0.1;
            q(99) = 0.1;
            q(97) = 0.1;
    end
else
    q = ones(1,length(alpha1));
    p = ones(1,length(alpha1));
    
end


pos_night = [40, 150, 40];
pos_day = [-60, -250, -60];

for i_reg = 1:3

label = {'','over_est','under_est'};

eval(['id_reg = id_' label{i_reg} ';'])

% obs:
sa_d_m = nanmean(day_sa_m(:,id_reg),2);
sa_n_m = nanmean(night_sa_m(:,id_reg),2);

%x_max_ = [5.5,17,4.5];


alpha1_reg = alpha1(id_reg);  % epi
alpha2_reg = alpha2(id_reg) .* q(id_reg)'; % shallow migrant
alpha5_reg = alpha2(id_reg) .* (1-q(id_reg)'); % deep migrant
alpha3_reg = alpha3(id_reg) .* p(id_reg)'; % shallow resident 
alpha4_reg = alpha3(id_reg) .* (1-p(id_reg)'); % deep resident


alpha_com = [];
for i=1:ng
    eval(['alpha' num2str(i) '_reg_ = repmat(alpha' num2str(i) '_reg,1,100);'])
    eval(['alpha_com = cat(3,alpha_com,alpha' num2str(i) '_reg_);'])
end
alpha_com = permute(alpha_com,[3,2,1]);


for l=1:length(len)
      d = importdata(['sigma_' num2str(l) '.txt']);
      phi = d.data; phi = repmat(phi,1,dn,ng,nx); phi = permute(phi,[2 3 1 4]);
      toto(:,:,:,:,l) = phi * ((len(l)/len(1))^(4*beta));
end



%load([pathway 'size_explicit/profile_norm_by_com_' project '_' prj '_now' suffix{1}]) %,'profile_norm'
load(['profile_by_com_' prj])
%%
load('depth')
id_depth = find(depth<=1000);
depth = -depth(id_depth);
if(size(depth,1)>1)
    depth = depth';
end
%%


profile_norm_ = squeeze(nansum(profile_by_com .* toto,5));
profile_norm_sum = squeeze(nansum(profile_norm_,3)); profile_norm_sum = repmat(profile_norm_sum,1,1,1,nz); profile_norm_sum = permute(profile_norm_sum,[1 2 4 3]);
profile_norm = profile_norm_./profile_norm_sum;


day_mod_com_m = squeeze(profile_norm(1,:,:,:)); %*100
night_mod_com_m = squeeze(profile_norm(2,:,:,:)); %*100

day_mod_com_m_ = day_mod_com_m(:,:,id_reg) .* alpha_com .* acoustic_day(:,:,id_reg);
night_mod_com_m_ = night_mod_com_m(:,:,id_reg) .* alpha_com .* acoustic_night(:,:,i_reg) ;

% Reverse order of communities 4 and 5
v = squeeze(day_mod_com_m_(5,:,:));
day_mod_com_m_(5,:,:) =  squeeze(day_mod_com_m_(4,:,:));
day_mod_com_m_(4,:,:) = v;

v = squeeze(night_mod_com_m_(5,:,:));
night_mod_com_m_(5,:,:) =  squeeze(night_mod_com_m_(4,:,:));
night_mod_com_m_(4,:,:) = v;


% mean (spatial)
day_mod_com_reg = nanmean(day_mod_com_m_,3); 
night_mod_com_reg = nanmean(night_mod_com_m_,3);
% quartiles
day_mod_com_reg_q1 = quantile(day_mod_com_m_,0.25,3);
day_mod_com_reg_q3 = quantile(day_mod_com_m_,0.75,3);
night_mod_com_reg_q1 = quantile(night_mod_com_m_,0.25,3);
night_mod_com_reg_q3 = quantile(night_mod_com_m_,0.75,3);

x_fill_d = [day_mod_com_reg_q1,day_mod_com_reg_q3(:,end:-1:1)];
x_fill_n = [night_mod_com_reg_q1,night_mod_com_reg_q3(:,end:-1:1)];

tot_day = nansum(day_mod_com_reg);
tot_night = nansum(night_mod_com_reg);

x_max = max(nanmax(tot_day),nanmax(tot_night));


eval(['ax = ax' num2str(i_reg) ';'])

fill(ax, [-x_max, -x_max 0 0],[max(depth) min(depth) min(depth) max(depth)],'w','facealpha',0.2)
hold(ax,"on") ;
fill(ax, [0 0 x_max x_max],[max(depth) min(depth) min(depth) max(depth)],[0.5 0.5 0.5],'facealpha',0.15)

% 

% daytime


xdata = [-tot_day,zeros(1,size(day_mod_com_reg,2))];

hold(ax,"on") ;

for com=1:5
    xdata = [-day_mod_com_reg(com,:),zeros(1,size(day_mod_com_reg,2))];
    ydata = [depth , depth(end:-1:1)];
    
    eval(['col = col' num2str(com) ';'])
    f(com) = fill(ax,xdata,ydata,col,'EdgeColor',col,'LineWidth',2,'FaceAlpha',0.4,'DisplayName',legend_lab{5-com+1});

    
end

hLines(1) = plot(ax,-tot_day,depth,'--k','linewidth',1,'DisplayName','Total');

% Night


xdata = [tot_night,zeros(1,size(day_mod_com_reg,2))];



hold(ax,"on") ;

for com=1:5
    xdata = [night_mod_com_reg(com,:),zeros(1,size(day_mod_com_reg,2))];
    ydata = [depth , depth(end:-1:1)];
    
    eval(['col = col' num2str(com) ';'])
    fill(ax,xdata,ydata,col,'EdgeColor',col,'LineWidth',2,'FaceAlpha',0.4)

    hold(ax,"on") ;
end

hLines(2) = plot(ax,tot_night,depth,'--k','linewidth',1);


pas = 2000;
if(i_reg==2)
    pas=10000;
end

text(ax,pos_day(i_reg),-850,'Day','fontsize',14)
text(ax,pos_night(i_reg),-850,'Night','fontsize',14)

x_tick = get(ax,'xtick');
%x_tick = min(x_tick):pas:max(x_tick);
x_tick_ = abs(x_tick);
x_tick_label = strsplit(num2str(x_tick_));
set(ax,'xtick',x_tick,'xticklabel',x_tick_label,'fontsize',12)

xlabel(ax,"Acoustic backscattering intensity","FontSize",14)

title(ax,titre{i_reg},'Color',couleur{i_reg}, 'FontSize',18)

if(i_reg == 3)
    
    lgd = legend([f(1:5),hLines(1)]);
    lgd.Location = "northwest";
    lgd.Box="off" ;
    lgd.FontSize = 11;
end

if(i_reg == 1)
    ylabel(ax,"Depth (m)","FontSize",14)
end

xlim(ax,[-x_max x_max])
ylim(ax,[-900 depth(1)])



end

set(ax1,'YTick',-900:100:0,'YTickLabel',{})
set(ax2,'YTick',-900:100:0,'YTickLabel',{'','-800','','-600','','-400','','-200','','0'})
set(ax3,'YTick',-900:100:0,'YTickLabel',{'','-800','','-600','','-400','','-200','','0'})


sgtitle('Contribution of the 5 communities to the total acoustic backscatter intensity ','FontSize',20)



