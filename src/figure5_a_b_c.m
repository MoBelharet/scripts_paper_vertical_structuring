%clearvars;

addpath("data/")

%%
fig = figure('Renderer', 'painters', 'Position', [10 10 1000 400]);



ax1 = subplot_tight(1,3,1,[0.2, 0.05]);

set(ax1,'box' , 'on');

ax2 = subplot_tight(1,3,2,[0.2, 0.05]);
set(ax2,'box' , 'on');

ax3 = subplot_tight(1,3,3,[0.2, 0.05]);
set(ax3,'box' , 'on');

%%

d = importdata('data/cluster_list_selected.txt');
clstr = d.data;

load('data/day_sa_m')
load('data/night_sa_m')

day_sa = day_sa_m(:,clstr);
night_sa = night_sa_m(:,clstr);



run('data_for_boxplot_paper.m');

%%

load('id_')
load('id_over_est')
load('id_under_est')

label = {'','over_est','under_est'};
titre = {'(a) GP profiles','(b) SP profiles','(c) DP profiles'};

leg = {'Simple' , 'Full', '', 'Simple' , 'Full';...
'Model' , 'Model', '', 'Model' , 'Model'};


%leg = {'Simple Model','Full Model','','Simple Model','Full Model'};

col_id = {[75 0 130]/225,[0 86 27]/225,[1 0 0]};

val = [0.13, 0.13, 0.09];

for i_reg = 1:3

eval(['id_reg = id_' label{i_reg} ';'])



%%

v1 = r_square_d(2,id_reg)';
v2 = r_square_d(1,id_reg)';
v3 = nan(length(id_reg),1);
v4 = r_square_n(2,id_reg)';
v5 = r_square_n(1,id_reg)';
%------ stats ---------------------
for i=1:5
    

    if(i~=3)
        eval(['v' num2str(i) '_min = quantile(v' num2str(i) ',0.1);'])
        eval(['v' num2str(i) '_max = quantile(v' num2str(i) ',0.9);'])

        eval(['v' num2str(i) '_median = quantile(v' num2str(i) ',0.5);'])
    else
        eval(['v' num2str(i) '_min = nan;'])
        eval(['v' num2str(i) '_max = nan;'])

        eval(['v' num2str(i) '_median = nan;'])

    end

end



%%


    eval(['ax = ax' num2str(i_reg) ';'])

    fill(ax,[3 5.5 5.5 3],[0 0 1 1],[0.5 0.5 0.5],'facealpha',0.1)

    hold(ax,'on')

    fill(ax,[0.5 3 3 0.5],[0 0 1 1],'w')
    set(ax,'xlim',[0.5 5.5])


    %L_bxp = boxplot(ax,[v1,v2,v3, v4,v5],'symbol','','labels',leg);
    L_bxp = boxplot(ax,[v1,v2,v3,v4,v5],'Symbol','');
   
    tickCell = {'XTickLabel',{}};
    set(ax,tickCell{:})

    for ii = 1:length(leg)
        h(ii)=text(ax, ax.XTick(ii), ax.YLim(1)-val(i_reg), sprintf('%s\n%s\n', leg{:,ii}), ...
            'horizontalalignment', 'center', 'verticalalignment', 'middle');
        set(h(ii),'Rotation',60); %rotating the x label
    end


    lines = ax.Children;

  uw = findobj(lines, 'tag', 'Upper Whisker');           % get handle to "Upper Whisker" line
  uav = findobj(lines, 'tag', 'Upper Adjacent Value');   %get handle to "Upper Adjacent Value" line
  lw = findobj(lines, 'tag', 'Lower Whisker');           % get handle to "Lower Whisker" line
  lav = findobj(lines, 'tag', 'Lower Adjacent Value');   %get handle to "Lower Adjacent Value" line

  box = findobj(lines, 'tag','Box'); %get handle to box
  med = findobj(lines, 'tag','Median'); %get handle to median

    for i=1:5
        j = 5-i+1;
        if(i~=3)
            eval(['uw(' num2str(i) ').YData(1,2) = v' num2str(j) '_max;'])
            eval(['uav(' num2str(i) ').YData(:) = v' num2str(j) '_max;'])
            eval(['lw(' num2str(i) ').YData(1,1) = v' num2str(j) '_min;'])
            eval(['lav(' num2str(i) ').YData(:) = v' num2str(j) '_min;'])

            % Color
            eval(['uw(' num2str(i) ').Color = col_id{i_reg};'])
            eval(['uav(' num2str(i) ').Color = col_id{i_reg};'])
            eval(['lw(' num2str(i) ').Color = col_id{i_reg};'])
            eval(['lav(' num2str(i) ').Color = col_id{i_reg};'])
            
            eval(['box(' num2str(i) ').Color = col_id{i_reg};'])
            eval(['med(' num2str(i) ').Color = col_id{i_reg};'])

            % Line style
            
            eval(['uw(' num2str(i) ').LineStyle =' ''':''' ';'])
            eval(['uw(' num2str(i) '). LineWidth = 1;' ]) 
            eval(['lw(' num2str(i) ').LineStyle =' ''':''' ';'])
            eval(['lw(' num2str(i) '). LineWidth = 1;' ])

            % Line width
            eval(['box(' num2str(i) ').LineWidth = 1.5;'])
            eval(['med(' num2str(i) ').LineWidth = 1.5;'])
        end
        

    end

  
    title(ax,titre{i_reg},'FontSize',12,'FontWeight','bold')

    text(ax,1.25,0.08,'Day','FontSize',12)
    text(ax,4.1,0.08,'Night','FontSize',12)


end

sgtitle('Model improvement per group','FontSize',13,'FontWeight','bold')

ylabel(ax1,'R²','fontsize',14)
ylim(ax1,[0 1])


set(ax1,'YTick',0:0.1:1)
tickCell = {'YTickLabel',{}};
set(ax1,tickCell{:})

set(ax2,'YTick',0:0.1:1,'YTickLabel',{'0','','0.2','','0.4','','0.6','','0.8','','1.0'})
set(ax2,'xlim',[0.5 5.5])
ylim(ax2,[0 1])


set(ax3,'YTick',0:0.1:1, 'YTickLabel',{'0','','0.2','','0.4','','0.6','','0.8','','1.0'})
set(ax3,'xlim',[0.5 5.5])
ylim(ax3,[0 1])

%%

%print(fig,'figures/fig6_abc','-dpng')