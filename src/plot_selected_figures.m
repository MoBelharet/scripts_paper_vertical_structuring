

%% ------------------------------------------------
if(plot_fig2)
    plot_fig2_a_b = 1;
    plot_fig2_c_d_e = 1;
else
    plot_fig2_a_b = 0;
    plot_fig2_c_d_e = 0;
end

if(plot_fig3)
    plot_fig3_a_b = 1;
    plot_fig3_c = 1;
    plot_fig3_d = 1;
else
    plot_fig3_a_b = 0;
    plot_fig3_c = 0;
    plot_fig3_d = 0;
end

if(plot_fig4)
    plot_fig4_a_b = 1;
    plot_fig4_c_d_e = 1;
else
    plot_fig4_a_b = 0;
    plot_fig4_c_d_e = 0;
end

if(plot_fig5)
    plot_fig5_a_b_c = 1;
    plot_fig5_d_e_f_g_h_i = 1;
else
    plot_fig5_a_b_c = 0;
    plot_fig5_d_e_f_g_h_i = 0;
end

%%
%***************************************************
%  Figure1
%***************************************************
if(plot_fig1)
    clearvars -except plot_fig* add_simple_model;
    addpath data/
    addpath functions/
    addpath src/
    prj = 'three_groups_without_oxy';
    figure1;
end
%%
%***************************************************
%  Figure2
%***************************************************
if(plot_fig2_a_b)
    clearvars -except plot_fig* add_simple_model;
    addpath data/
    addpath functions/
    addpath src/
    prj = 'three_groups_without_oxy';
    scatter_plots;
end


if(plot_fig2_c_d_e)
    clearvars -except plot_fig* add_simple_model;
    addpath data/
    addpath functions/
    addpath src/

    add_simple_model = 0;
    plot_fig3_a_b_ = plot_fig3_a_b;
    plot_fig3_c_ = plot_fig3_c;
    plot_fig3_a_b = 0;
    plot_fig3_c = 0;
    prj = 'three_groups_without_oxy';
    figure2_c_d_e_figure3_a_b;
    plot_fig3_a_b = plot_fig3_a_b_;
    plot_fig3_c = plot_fig3_c_;
end
%****************************************************

%%
%***************************************************
% Figure3
%***************************************************

if(plot_fig3_a_b)
    clearvars -except plot_fig* add_simple_model ;
    addpath data/
    addpath functions/
    addpath src/

    plot_fig2_c_d_e_ = plot_fig2_c_d_e;
    plot_fig2_c_d_e = 0;
    plot_fig3_c_ = plot_fig3_c;
   
    prj = 'three_groups_without_oxy';
    plot_fig3_c = 0;
    add_simple_model=0;
    figure2_c_d_e_figure3_a_b;

    plot_fig2_c_d_e = plot_fig2_c_d_e_;
    plot_fig3_c = plot_fig3_c_;
    
end


if(plot_fig3_c)
    clearvars -except plot_fig* add_simple_model ;
    addpath data/
    addpath functions/
    addpath src/

    plot_fig3_a_b_ = plot_fig3_a_b;
    
    plot_fig2_c_d_e_ = plot_fig2_c_d_e;

    plot_fig3_a_b = 0;
    

    prj = 'three_groups_with_oxy';
    add_simple_model=1;
    figure2_c_d_e_figure3_a_b;

    plot_fig3_a_b = plot_fig3_a_b_;
    plot_fig2_c_d_e = plot_fig2_c_d_e_;
end


if(plot_fig3_d)
    clearvars -except plot_fig* add_simple_model;
    addpath data/
    addpath functions/
    addpath src/


    prj = 'three_groups_with_oxy';
    scatter_plots;
end
%***************************************************

%%
%***************************************************
% Figure4, 
%***************************************************
if(plot_fig4_a_b)
    clearvars -except plot_fig* add_simple_model;
    addpath data/
    addpath functions/
    addpath src/

    prj ='five_groups_with_oxy';%
    scatter_plots;
end


if(plot_fig4_c_d_e)
    clearvars -except plot_fig* add_simple_model;
    addpath data/
    addpath functions/
    addpath src/

    plot_fig3_a_b_ = plot_fig3_a_b;
    plot_fig3_c_ = plot_fig3_c;

    plot_fig2_c_d_e_ = plot_fig2_c_d_e;

    plot_fig3_a_b = 0;
    plot_fig3_c = 0;
    plot_fig2_c_d_e = 1;

    prj = 'five_groups_with_oxy';
    add_simple_model = 1;
    figure2_c_d_e_figure3_a_b;
    %close all
    
    plot_fig3_a_b = plot_fig3_a_b_;
    plot_fig3_c = plot_fig3_c_;

    plot_fig2_c_d_e = plot_fig2_c_d_e_;
end
%****************************************************
%%

%***************************************************
% Figure5,
%***************************************************

if(plot_fig5_a_b_c)
clearvars -except plot_fig* add_simple_model;
    addpath data/
    addpath functions/
    addpath src/
    figure5_a_b_c;

end

if(plot_fig5_d_e_f_g_h_i)
clearvars -except plot_fig* add_simple_model ;
    addpath data/
    addpath functions/
    addpath src/
    figure5_d_e_f_g_h_i;

end

%%
%***************************************************
% Figure6,
%***************************************************
if(plot_fig6)
    clearvars -except plot_fig* add_simple_model;
    addpath data/
    addpath functions/
    addpath src/

    
    figure6

    
end
