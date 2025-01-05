
%clearvars;
addpath("data/")
addpath('data/othercolor/')

%%

id_oc = sort([15 30 55 71 82 100]);
id_oc_mean = mean([[0,id_oc]; [id_oc, 122]]);
 
leg = {'North', 'South', 'South', 'South', 'Nwest.','Neast.','North' ;...
    'Atlantic', 'Atlantic', 'Indian', 'Pacific','Pacific.','Pacific','Atlantic'};

%%
fig = figure('Renderer', 'painters', 'Position', [10 10 1000 800]); 


for i=1:6

    eval(['ax' num2str(i) '= subplot_tight(3,2,' num2str(i) ',[0.08, 0.05]);']) %[0.08, 0.05]

    eval(['set(ax' num2str(i) ',' '''box''' ',' '''on''' ');'])

end




%%
d = importdata('cluster_list_selected.txt');
clstr_list = d.data;

%lat_m = lat_m (clstr_list);
%lon_m = lon_m (clstr_list);

load('day_sa_m.mat')
load('night_sa_m.mat')

load('depth.mat')
depth = depth(1:100);
if(size(depth,2)>1)
    depth = depth';
end

[clstr_list_sorted, id_sort] = sort(clstr_list,'ascend');
id_depth = 1:90;

%%

prj = 'three_groups_without_oxy';
% Model
load(['profile_total_MALASPINA_' prj ])

profile_total = quantile(profile_total,0.5,4);
profile_total_d = squeeze(profile_total(1,:,:,:));
profile_total_n = squeeze(profile_total(2,:,:,:));

% Obs
day_sa_m = day_sa_m(id_depth,:)./nansum(day_sa_m(id_depth,:));
night_sa_m = night_sa_m(id_depth,:)./nansum(night_sa_m(id_depth,:));
sa_d = day_sa_m(:,clstr_list);
sa_n = night_sa_m(:,clstr_list);
sa_d = sa_d(:,id_sort);
sa_n = sa_n(:,id_sort);
sa_d = sa_d./repmat(nansum(sa_d),length(id_depth),1);
sa_n = sa_n./repmat(nansum(sa_n),length(id_depth),1);


DEPTH  = repmat(depth(id_depth),1,size(profile_total_d,2));
X = repmat(1:size(profile_total_d,2),length(id_depth),1);


%%

idx = [1,20:20:100, 120];

colorscheme = othercolor('Paired8',50);


colormap(colorscheme)

%%

% Day
pcolor(ax1,X(id_depth,:),-DEPTH(id_depth,:),profile_total_d(id_depth,:)); shading(ax1,"interp") ;
hold(ax1,"on") 
caxis(ax1,[0 0.08]);


for i=1:length(id_oc)
    plot(ax1, [id_oc(i) id_oc(i)],[-900 1],'--k','LineWidth',0.5)
end
%set(ax1,'fontsize',9,'fontweight','bold')
set(ax1,'XTick',[],'YTick',[] )
ylabel(ax1,'Depth (m)','FontSize',12)

pos = get(ax1, 'position');
dim = cellfun(@(x) x.*[5.95 0.9 0.5 0.48], mat2cell(pos,1), 'uni',0);
annotation(fig, 'textbox', dim{1}, 'String',"(d) Day : Simple model " , 'FitBoxToText','on','BackgroundColor','w');

for ii = 1:length(leg)
      h(ii)=text(ax1, id_oc_mean(ii), ax1.YLim(1)-150, sprintf('%s\n%s\n', leg{:,ii}), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle');
      set(h(ii),'Rotation',60); %rotating the x label
end

%annotation(ax1, 'textbox', [0.5, 0.2, 0.1, 0.1], 'String', "(d) Day : Simple model ")

% Night

pcolor(ax2,X(id_depth,:),-DEPTH(id_depth,:),profile_total_n(id_depth,:)); shading(ax2,"interp") ;
hold(ax2,"on") 
caxis(ax2,[0 0.08]);



for i=1:length(id_oc)
    plot(ax2, [id_oc(i) id_oc(i)],[-900 1],'--k','LineWidth',0.5)
end
%set(ax1,'fontsize',9,'fontweight','bold')
set(ax2,'XTick',[])

pos = get(ax2, 'position');
dim = cellfun(@(x) x.*[1.4542 0.9 0.5 0.48], mat2cell(pos,1), 'uni',0);
annotation(fig, 'textbox', dim{1}, 'String',"(e) Night : Simple model " , 'FitBoxToText','on','BackgroundColor','w');

for ii = 1:length(leg)
      h(ii)=text(ax2, id_oc_mean(ii), ax2.YLim(1)-150, sprintf('%s\n%s\n', leg{:,ii}), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle');
      set(h(ii),'Rotation',60); %rotating the x label
end

%%
prj = 'five_groups_with_oxy';

load(['profile_total_MALASPINA_' prj ])

profile_total = quantile(profile_total,0.5,4);
profile_total_d = squeeze(profile_total(1,:,:,:));
profile_total_n = squeeze(profile_total(2,:,:,:));

pcolor(ax3,X(id_depth,:),-DEPTH(id_depth,:),profile_total_d(id_depth,:)); shading(ax3,"interp") ;
hold(ax3,"on") 
caxis(ax3,[0 0.08]);


for i=1:length(id_oc)
    plot(ax3, [id_oc(i) id_oc(i)],[-900 1],'--k','LineWidth',0.5)
end
%set(ax3,'fontsize',9,'fontweight','bold')
set(ax3,'XTick',[],'YTick',[] )
ylabel(ax3,'Depth (m)','FontSize',12)

pos = get(ax3, 'position');
dim = cellfun(@(x) x.*[6.45 0.85 0.5 0.42], mat2cell(pos,1), 'uni',0);
annotation(fig, 'textbox', dim{1}, 'String',"(f) Day : Full model " , 'FitBoxToText','on','BackgroundColor','w');

for ii = 1:length(leg)
      h(ii)=text(ax3, id_oc_mean(ii), ax3.YLim(1)-150, sprintf('%s\n%s\n', leg{:,ii}), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle');
      set(h(ii),'Rotation',60); %rotating the x label
end


% Night

pcolor(ax4,X(id_depth,:),-DEPTH(id_depth,:),profile_total_n(id_depth,:)); shading(ax4,"interp") ;
hold(ax4,"on") 
caxis(ax4,[0 0.08]);


for i=1:length(id_oc)
    plot(ax4, [id_oc(i) id_oc(i)],[-900 1],'--k','LineWidth',0.5)
end

set(ax4,'XTick',[])
%ylabel(ax4,'Depth (m)','FontSize',12)

pos = get(ax4, 'position');
dim = cellfun(@(x) x.*[1.496 0.85 0.5 0.42], mat2cell(pos,1), 'uni',0);
annotation(fig, 'textbox', dim{1}, 'String',"(g) Night : Full model " , 'FitBoxToText','on','BackgroundColor','w');

for ii = 1:length(leg)
      h(ii)=text(ax4, id_oc_mean(ii), ax4.YLim(1)-150, sprintf('%s\n%s\n', leg{:,ii}), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle');
      set(h(ii),'Rotation',60); %rotating the x label
end

%%

% Day
pcolor(ax5,X(id_depth,:),-DEPTH(id_depth,:),sa_d(id_depth,:)); shading(ax5,"interp");
hold(ax5,"on") 
caxis(ax5,[0 0.08])

for i=1:length(id_oc)
    plot(ax5,[id_oc(i) id_oc(i)],[-900 1],'--k','LineWidth',0.5)
end
set(ax5,'XTick',[],'YTick',[] )
ylabel(ax5,'Depth (m)','FontSize',12)


pos = get(ax5, 'position');
dim = cellfun(@(x) x.*[7.24 0 0.5 0.52], mat2cell(pos,1), 'uni',0);
annotation(fig, 'textbox', dim{1}, 'String',"(h) Day : Obs " , 'FitBoxToText','on','BackgroundColor','w');

for ii = 1:length(leg)
      h(ii)=text(ax5, id_oc_mean(ii), ax5.YLim(1)-150, sprintf('%s\n%s\n', leg{:,ii}), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle');
      set(h(ii),'Rotation',60); %rotating the x label
end


% Night

pcolor(ax6,X(id_depth,:),-DEPTH(id_depth,:),sa_n(id_depth,:)); shading(ax6,"interp");
hold(ax6,"on") 
caxis(ax6,[0 0.08])

for i=1:length(id_oc)
    plot(ax6,[id_oc(i) id_oc(i)],[-900 1],'--k','LineWidth',0.5)
end
set(ax6,'XTick',[] )
%ylabel(ax5,'Depth (m)','FontSize',12)


pos = get(ax6, 'position');
dim = cellfun(@(x) x.*[1.585 0 0.5 0.52], mat2cell(pos,1), 'uni',0);
annotation(fig, 'textbox', dim{1}, 'String',"(i) Night : Obs " , 'FitBoxToText','on','BackgroundColor','w');


for ii = 1:length(leg)
      h(ii)=text(ax6, id_oc_mean(ii), ax6.YLim(1)-150, sprintf('%s\n%s\n', leg{:,ii}), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle');
      set(h(ii),'Rotation',60); %rotating the x label
end


%%
sgtitle('Profiles around Malaspina cruise track','FontSize',15,'FontWeight','bold')
%print(fig,'figures/fig6_defghi','-dpng')
