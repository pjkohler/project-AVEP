fig_path = '/Users/kohler/Google Drive/WRITING/Articles/2018_AVEP/figures';

load(sprintf('%s/avep_sums.mat',fig_path));
n_groups = size(full_sum, 1);
n_conds = size(full_sum, 2);

f_size = 12;
l_width = 1;
gca_opts = {'tickdir','out','ticklength',[0.0200,0.0200],'box','off','fontsize',f_size,'fontname','Helvetica','linewidth',l_width,'clipping','off'};

close all

fig_h = figure;
set(fig_h,'units','centimeters');
fig_pos = get(gcf,'pos');
fig_pos(3) = 15; fig_pos(4) = 150;
set(gcf,'pos',fig_pos);
sub_groups = {'TD','ASD','ADHD'};
cond_names = {'vis sweep','aud sweep','vis sweep + aud','aud sweep + vis'};


for c = 1:n_conds
    for g = 1:n_groups
        figure(fig_h);
        s_h = subplot(n_conds, length(subject_group), g+n_groups*(c-1));
        hold on;
        imagesc(full_sum{g,c}, [.6, 1]);
        min_val(g,c) = min(full_sum{g,c}(:));
        median_val(g,c) = median(full_sum{g,c}(:));
        xlim([0.5,10.5])
        ylim([0.5,128.5])
        set(gca, gca_opts{:}, 'xtick', 1:10, 'xticklabels', {'','2','','4','','6','','8','','10'}, 'ytick', [1,10:10:120, 128]);
        if c == 1
            title(['\it{',sub_groups{g},'}'], 'fontsize',f_size,'fontname','Helvetica', 'fontweight','normal', 'horizontalalignment','center');
        else
        end
        if g == 1
            if c == n_conds
                xlabel('\it{bins}', 'fontsize',f_size,'fontname','Helvetica')
            else
            end
            ylabel(['\it{',cond_names{c},'}'], 'fontsize',f_size,'fontname','Helvetica')
        else
        end
        if c == n_conds && g == 2
            colorbar( 'southoutside', 'fontsize',f_size,'fontname','Helvetica','linewidth', l_width, 'xtick', .6:.2:1)
            new_pos = get(s_h, 'position');
            fig_pos(1) = new_pos(1); 
        else
            fig_pos = get(s_h, 'position');
            fig_pos(2) = fig_pos(2)+fig_pos(4)*.2*c;
        end
        set(s_h, 'position', fig_pos);
        hold off
    end
end
export_fig(sprintf('%s/avep_sums.png',fig_path), '-png', '-r600', '-transparent', fig_h);
