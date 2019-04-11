codeFolder = '/Users/kohler/code';
rcaCodePath = sprintf('%s/git/rcaBase', codeFolder);
addpath(genpath(rcaCodePath));
addpath(genpath(sprintf('%s/git/mrC', codeFolder)));
addpath(genpath(sprintf('%s/git/schlegel/matlab_lib/figure', codeFolder)));
addpath(genpath(sprintf('%s/git/schlegel/matlab_lib/misc', codeFolder)));

addpath(genpath(sprintf('%s/git/export_fig', codeFolder)));
addpath(genpath(sprintf('%s/git/mes_toolbox',codeFolder)));

setenv('DYLD_LIBRARY_PATH','')
close all;
clear all;

printFigures=true; % set to true if you want to automatically print the figures and save them into the high-level directory where the data are
run_rca = true;
trial_error = false;
plot_medians = false;
plot_binmean = false;

bins_to_use =1:10; % indices of bins to include in analysis (the values must be present in the bin column of all DFT/RLS exports)
trials_to_use = []; %1:10; % subset of trials to use for analysis (if set to false or empty, all trials will be used)
nReg=7; % RCA regularization constant (7-9 are typical values, but see within-trial eigenvalue plot in rca output)
nComp=5; % number of RCs that you want to look at (3-5 are good values, but see across-trial eigenvalue plot in rca output)
chanToCompare = 75; % channel to use for a performance evaluation, can be []
data_type = 'DFT'; % can also be 'DFT' if you have DFT exports
rcPlotStyle = 'matchMaxSignsToRc1'; % not req'd. see 'help rcaRun', can be: 'matchMaxSignsToRc1' (default) or 'orig'
forceSourceData = false;
doNR = false(4,6,5);  % 10 freqs, 5 RCs, 5 conditions
%doNR(:,:,:) = true; % do fitting for all RC, first harmonic, all conditions

%% GET SUBJECT FOLDERS
top_path = '/Volumes/Denali_DATA1/kohler/EEG_EXP/DATA/AVEP';
subject_group{1} = sprintf('%s/AVEP_ASD',top_path); % ASD young subjects
subject_group{2} = sprintf('%s/AVEP_TYP',top_path); % Typical young subjects
subject_group{3} = sprintf('%s/AVEP_ADHD',top_path);% ADHD young participants
sub_count = 0;
for d=1:length(subject_group)
    folder_names = subfolders(sprintf('%s/SP*',subject_group{d}),1);
    wobd_idx = cell2mat(cellfun(@(x) ~isempty(strfind(x,'wobd')),folder_names,'uni',false));
    folder_names = folder_names(~wobd_idx);
    for f = 1:length(folder_names)
        temp_names = subfolders(sprintf('%s/*Exp_TEXT_*',folder_names{f}),1);
        if temp_names{1} == 0
            msg = sprintf('missing data for subject %s',folder_names{f}); 
            warning(msg);
        else
            sub_count = sub_count + 1;
            ready_names{sub_count} = temp_names{1};
        end
    end
end
asd_idx = cell2mat(cellfun(@(x) ~isempty(strfind(upper(x),'ASD')),ready_names,'uni',false));
typ_idx = cell2mat(cellfun(@(x) ~isempty(strfind(upper(x),'TYP')),ready_names,'uni',false));
adhd_idx = cell2mat(cellfun(@(x) ~isempty(strfind(upper(x),'ADHD')),ready_names,'uni',false));

%% GET THE RELEVANT DATA, AND RUN RCA
% we are running RCA twice, on all subjects
% first, with 1F1 and 2F1, using the visual condition (1)
% then, with 1F2 and 2F2, using the auditory condition (2)
% we want to apply the components from each of these RCS
% to the audio-visual conditions (3,4) 

cond_names = {'vis sweep','aud sweep','vis sweep + aud','aud sweep + vis'};
all_freqs = [1,2,5,6];
all_conds = [1,2,3,4];
mat_path = sprintf('%s/analysis/rca_avep.mat',top_path);
errorType = 'SEM';
keepConditions = true;

if ~run_rca % if users says analysis is not to be launched
    load(mat_path);
else
    for q = 1:2
        if q == 1
            freq_idx = [1,2]; 
            conds_to_use = all_conds([1,3]);
            force_source_data = false; % true
            rca_figpath = sprintf('%s/analysis/rca_vis.pdf',top_path);
        else
            freq_idx = [3,4]; 
            conds_to_use = all_conds([2,4]);
            force_source_data = false;
            rca_figpath = sprintf('%s/analysis/rca_aud.pdf',top_path);
        end
        rca_all(freq_idx) = rcaSweep(ready_names,bins_to_use,all_freqs(freq_idx),conds_to_use,trials_to_use,nReg,nComp,data_type,chanToCompare,[],rcPlotStyle,force_source_data);
        export_fig(rca_figpath,'-pdf','-transparent','-opengl',gcf);
    end
    save(mat_path,'rca_all'); 
end

%% NOW PROJECT ALL DATA THROUGH THE TWO COMPONENTS, SPLITTING BY GROUP
warning('off','all')
for g = 1:length(subject_group)
    if g == 1
        group_names = ready_names(typ_idx);
    elseif g == 2
        group_names = ready_names(asd_idx);
    else
        group_names = ready_names(adhd_idx);
    end  
    for s = 1:length(group_names)
        [signal_data(:,s),noise_data1(:,s),noise_data2(:,s),sub_freq_idx{s},sub_bin_idx{s},sub_freq_labels{s},sub_bin_labels{s}] = selectDataForTraining(sprintf('%s/sourceData_%s.mat',group_names{s},data_type),bins_to_use,all_freqs,[],trials_to_use);
    end
    % CHECK FREQUENCY AND BIN INDICES FOR CONSISTENCY ACROSS SUBS
    n_subs(g) = size(signal_data,2);
    for s=1:n_subs(g)
        for c = 1:length(all_conds)
            if sum(abs(sub_freq_idx{s}{c}-sub_freq_idx{1}{c}))~=0 && sum(abs(sub_bin_idx{s}{c}-sub_bin_idx{1}{c}))~=0
                error('Frequency and bin indices vary across subjects: check consistency of DFT/RLS exports\n.');
            end
        end
    end
    freq_labels = sub_freq_labels{1}{1};
    vis_labels = cellfun(@(x) num2str(str2num(x),'%0.2f'), sub_bin_labels{1}{1},'uni',false);
    aud_labels = cellfun(@(x) num2str(str2num(x),'%0.4f'), sub_bin_labels{1}{2},'uni',false);
    % GENERATE NEW RCA STRUCT WITH ALL GROUP-LEVEL DATA
    for q = 1:2
        for f = 1:length(all_freqs)
            if q == 1
                % 1 and 2 will contain rc comps based on vis
                projected(g,q).rca(f).W = rca_all(1).W;
                projected(g,q).rca(f).A = rca_all(1).A;
                projected(g,q).rca(f).settings = rca_all(1).settings;
            else
                % 3 and 4 will contain rc comps based on aud
                projected(g,q).rca(f).W = rca_all(3).W;
                projected(g,q).rca(f).A = rca_all(3).A;
                projected(g,q).rca(f).settings = rca_all(3).settings;
            end
    
            temp_rcasignal = rcaProject( signal_data, projected(g,q).rca(f).W );
            temp_rcanoise.lowerSideBand=rcaProject( noise_data1, projected(g,q).rca(f).W ); 
            temp_rcanoise.higherSideBand=rcaProject( noise_data2, projected(g,q).rca(f).W );

            % create a "component" of just one channel for performance evaluation if requested
            chanToCompare = 75;
            nChannels = size(signal_data{1},2);
            wComparison=zeros(nChannels,1); wComparison(chanToCompare)=1; 
            temp_comparison_signal = rcaProject(signal_data,wComparison); 
            temp_comparison_noise.lowerSideBand = rcaProject(noise_data1,wComparison); 
            temp_comparison_noise.higherSideBand = rcaProject(noise_data2,wComparison);
            
            % use freq_idx from first subject, they are all the same
            f_idx = repmat(sub_freq_idx{1}{1}==all_freqs(f),2,1); % repmat because the first half is real, second half is imag with same ordering
            
            projected(g,q).rca(f).data = cellfun(@(x) x(f_idx,:,:),temp_rcasignal,'uni',false);
            projected(g,q).rca(f).noiseData.lowerSideBand = cellfun(@(x) x(f_idx,:,:),temp_rcanoise.lowerSideBand,'uni',false);
            projected(g,q).rca(f).noiseData.higherSideBand = cellfun(@(x) x(f_idx,:,:),temp_rcanoise.higherSideBand,'uni',false);
            projected(g,q).rca(f).comparisonData = cellfun(@(x) x(f_idx,:,:),temp_comparison_signal,'uni',false);
            projected(g,q).rca(f).comparisonNoiseData.lowerSideBand = cellfun(@(x) x(f_idx,:,:),temp_comparison_noise.lowerSideBand,'uni',false);
            projected(g,q).rca(f).comparisonNoiseData.higherSideBand = cellfun(@(x) x(f_idx,:,:),temp_comparison_noise.higherSideBand,'uni',false);
            projected(g,q).rca(f).settings = orderfields(projected(g,q).rca(f).settings);
            temp_struct = aggregateData(projected(g,q).rca(f),keepConditions,errorType,trial_error,doNR);        
            % RC
            projected(g,q).rca(f).stats.Amp = squeeze(temp_struct.ampBins);
            projected(g,q).rca(f).stats.SubjectAmp = squeeze(temp_struct.subjectAmp);
            projected(g,q).rca(f).stats.SubjectAmpProjected = squeeze(temp_struct.subjectAmp_projected);
            projected(g,q).rca(f).stats.ErrLB = squeeze(temp_struct.ampErrBins(:,:,:,:,1));
            projected(g,q).rca(f).stats.ErrUB = squeeze(temp_struct.ampErrBins(:,:,:,:,2));
            projected(g,q).rca(f).stats.NoiseAmp = squeeze(temp_struct.ampNoiseBins);
            projected(g,q).rca(f).stats.SubjectNoiseAmp = squeeze(temp_struct.subjectAmpNoise);
            % Naka-Rushton
            projected(g,q).rca(f).stats.NR_Params = squeeze(temp_struct.NakaRushton.Params);
            projected(g,q).rca(f).stats.NR_R2 = squeeze(temp_struct.NakaRushton.R2);
            projected(g,q).rca(f).stats.NR_JKSE = squeeze(temp_struct.NakaRushton.JackKnife.SE);
            projected(g,q).rca(f).stats.NR_JKParams = squeeze(temp_struct.NakaRushton.JackKnife.Params);
            projected(g,q).rca(f).stats.hModel = temp_struct.NakaRushton.hModel;
            projected(g,q).rca(f).stats.tSqrdP = squeeze(temp_struct.tSqrdP);
            projected(g,q).rca(f).stats.tSqrdSig = squeeze(temp_struct.tSqrdSig);
            projected(g,q).rca(f).stats.tSqrdVal = squeeze(temp_struct.tSqrdVal);
            
            % make classication data
            temp_struct = aggregateData(projected(g,q).rca(f),keepConditions,errorType,true,false(4,6,5));
            temp_clfdata = squeeze(temp_struct.subjectAmp_projected);
            temp_clfdata = temp_clfdata(1:10,:,:,:); % don't include averages
            temp_targets = ones(size(temp_clfdata,3),1)*g;
            if g == 1
                clf_data{q,f} = temp_clfdata;
                clf_targets{q,f} = temp_targets;
            else
                clf_data{q,f} = cat(3,clf_data{q,f},temp_clfdata);
                clf_targets{q,f} = cat(1,clf_targets{q,f},temp_targets);
            end
            clear temp_*
        end       
    end
    clear signal* noise* sub_*;
end
warning('on','all')

%% DO CLASSIFICATION
% clf_q = 2;
% clf_c = 3;
% clf_rc = 1;
% clf_f = 1;
% 
% clf_t = clf_targets{clf_q,clf_f}; % nsample x nfeature
% clf_d = squeeze(clf_data{clf_q,clf_f}(:,clf_rc,:,clf_c))'; % nsample x nfeature
% clf_chunk = repmat(1:10,1,size(ready_names,2))'; % make trials separate chunks
% %clf_chunk(isnan(clf_d)) = 0;
% 
% res = MVPA.CrossValidation(clf_d,clf_t,clf_chunk,'partitioner',1,'target_balancer',10)

%% DO HOTELLING'S T2
for g = 1:length(subject_group)
    for q = 1:2
        for f = 1:length(all_freqs)
            [rca_data_real,rca_data_imag] = getRealImag(projected(g,q).rca(f).data);
             [comp_data_real,comp_data_imag] = getRealImag(projected(g,q).rca(f).comparisonData);
           
            % concatenate comparison and rca data
            rca_data_real = cellfun(@(x,y) cat(2,x,y), rca_data_real, comp_data_real, 'uni', false);
            rca_data_imag = cellfun(@(x,y) cat(2,x,y), rca_data_imag, comp_data_imag, 'uni', false);
            
            for r = 1:6
                x_data = cell2mat(permute(cellfun(@(x) squeeze(nanmean(x(:,r,:),3)), rca_data_real, 'uni', false),[3,2,1]));
                y_data = cell2mat(permute(cellfun(@(x) squeeze(nanmean(x(:,r,:),3)), rca_data_imag, 'uni', false),[3,2,1]));
                % compute bin average
                x_data = cat(1,x_data,nanmean(x_data,1));
                y_data = cat(1,y_data,nanmean(y_data,1));
                xy_data{g,q,f,r} = cat(2,permute(x_data,[2,4,1,3]),permute(y_data,[2,4,1,3])); % subjects x xy x bins x conditions
            end
        end
    end
end
% plotting variables
sub_groups = {'typ','asd','adhd'};
% colors
c_brewer = load('colorBrewer_new.mat');
group_colors = [c_brewer.rgb20(5,:); c_brewer.rgb20(7,:); c_brewer.rgb20(1,:)];
f_size = 12;
l_width = 1;
gcaOpts = {'tickdir','out','ticklength',[0.0400,0.0400],'box','off','fontsize',f_size,'fontname','Helvetica','linewidth',l_width,'clipping','off'};
%% do comparisons
use_projected = true;
t2_grp = [1,3]; t2_rc = 1; t2_frq = 1; t2_rctp = 1; t2_cnd = 1;
for b = 1:11
    tmp_strct = t2FC(xy_data{t2_grp(1),t2_rctp,t2_frq,t2_rc}(:,:,b,t2_cnd), xy_data{t2_grp(2),t2_rctp,t2_frq,t2_rc}(:,:,b,t2_cnd));
    between_tsqrd(b) = tmp_strct.tSqrdCritical;
    between_pval(b) = tmp_strct.pVal;
    between_mean1(b) = sqrt(nanmean(xy_data{t2_grp(1),t2_rctp,t2_frq,t2_rc}(:,1,b,t2_cnd))^2 + nanmean(xy_data{t2_grp(1),t2_rctp,t2_frq,t2_rc}(:,2,b,t2_cnd))^2 );
    between_mean2(b) = sqrt(nanmean(xy_data{t2_grp(2),t2_rctp,t2_frq,t2_rc}(:,1,b,t2_cnd))^2 + nanmean(xy_data{t2_grp(2),t2_rctp,t2_frq,t2_rc}(:,2,b,t2_cnd))^2 );
end
box_labels = cellfun(@(x) sprintf('%.1f',str2num(x)), vis_labels, 'uni', false);
box_labels = [repmat(box_labels,3,1)', {'ave'},  {'ave'},  {'ave'}];

if use_projected
    group1 = squeeze(projected(1,t2_rctp).rca(t2_frq).stats.SubjectAmpProjected(:,t2_rc,:,t2_cnd))';
    group2 = squeeze(projected(2,t2_rctp).rca(t2_frq).stats.SubjectAmpProjected(:,t2_rc,:,t2_cnd))';
    group3 = squeeze(projected(3,t2_rctp).rca(t2_frq).stats.SubjectAmpProjected(:,t2_rc,:,t2_cnd))';
    box_ylims = [-15,30];
else
    group1 = squeeze(sqrt(xy_data{1,t2_rctp,t2_frq,t2_rc}(:,1,:,t2_cnd).^2+xy_data{1,t2_rctp,t2_frq,t2_rc}(:,2,:,t2_cnd).^2));
    group2 = squeeze(sqrt(xy_data{2,t2_rctp,t2_frq,t2_rc}(:,1,:,t2_cnd).^2+xy_data{2,t2_rctp,t2_frq,t2_rc}(:,2,:,t2_cnd).^2));
    group3 = squeeze(sqrt(xy_data{3,t2_rctp,t2_frq,t2_rc}(:,1,:,t2_cnd).^2+xy_data{3,t2_rctp,t2_frq,t2_rc}(:,2,:,t2_cnd).^2));
    box_ylims = [0,35];
end
close all

%honest_plot(group2)
all_pos = [];
for z = 1:3
    if z == 1
        cur_data = group1;
    elseif z == 2
        cur_data = group2;
    else
        cur_data = group3;
    end
    bin_pos = (.5:.5:5)+(z-1)*5.5;
    boxplot(cur_data(:,1:10), 'positions', bin_pos , 'labels',1:10, 'width', .5, 'color', group_colors(z,:), 'symbol', 'o' )
    hold on
    boxplot(cur_data(:,11), 'positions', 16.5+z*.5, 'labels', 1, 'width', .5, 'color', group_colors(z,:), 'symbol', 'o' )
    all_pos = [all_pos, bin_pos];
end
set(gca,'ylim',box_ylims,'xlim',[0,19],'xtick', [all_pos, 17:.5:18], 'xticklabel', box_labels, 'clipping','off', gcaOpts{:},'ygrid', 'on')
set(gcf,'units','centimeters');
box_pos = get(gcf,'position');
box_pos(3) = 30; box_pos(4) = 10;
set(gcf,'position', box_pos);
text(max(get(gca,'xlim'))*.95,max(get(gca,'ylim')),sub_groups{1}, 'color', group_colors(1,:), 'fontsize',f_size*2,'fontname','Helvetica')
text(max(get(gca,'xlim'))*.95,max(get(gca,'ylim'))*.95,sub_groups{2}, 'color', group_colors(2,:), 'fontsize',f_size*2,'fontname','Helvetica')
text(max(get(gca,'xlim'))*.95,max(get(gca,'ylim'))*.9,sub_groups{3}, 'color', group_colors(3,:), 'fontsize',f_size*2,'fontname','Helvetica')

%% MAKE FIGURES
close all;

text_params = {'fontsize',f_size,'fontname','Helvetica','FontWeight','normal'};

for f = 1:4
    for q = 1:2
        fig_num = (f > 2)*2 + q;
        figure(fig_num);
        
        switch fig_num
            case 1
                y_min = 0;
                y_unit = 4;
                y_max = 20;
            case 2
                y_min = 0;
                y_unit = 4;
                y_max = 20;
            case 3
                y_min = 0;
                y_unit = 2;
                y_max = 8;
            otherwise
                y_min = 0;
                y_unit = 2;
                y_max = 8;
        end
        if ~mod(f,2)
            y_min = y_min/2;
            y_max = y_max/2;
            y_unit = y_unit/2;
        end
        str_array = {};
        for c = 1:4
            if mod(c,2)
                % x-values
                x_vals = cell2mat(cellfun(@(x) reallog(str2num(x)), vis_labels, 'uni',false));
                xtick_labels = [0.5,1,2,5,10,20,40,80];
                xticks = reallog(xtick_labels);
                bin_labels = vis_labels;
                x_str = 'contrast (%)';
            else
                x_vals = cell2mat(cellfun(@(x) reallog(str2num(x)), aud_labels, 'uni',false));
                x_vals = interp1([min(x_vals),max(x_vals)],[50,70],x_vals);
                xtick_labels = arrayfun(@(x) num2str(x,'%0.0f'), 50:5:70, 'uni',false);
                xticks = cellfun(@(x) str2num(x), xtick_labels);
                bin_labels = aud_labels;
                x_str = 'sound pressure (\it{dB})';
            end
            log_step = diff(x_vals(1:2)); % step size
            x_min = x_vals(1)-log_step*.5;
            if plot_binmean
                x_max = x_vals(end)+log_step*2.5; % add 2 steps
            else
                x_max = x_vals(end)+log_step*0.5;
            end
            extra_bins = arrayfun(@(x) x_vals(end)+x, [log_step,log_step*1.5,log_step*2]);
            s_h = subplot(2,5,c+~mod(f,2)*5);
            sub_pos = get(s_h,'position');
            sub_pos(4) = sub_pos(4) * 1.1;
            sub_pos = set(s_h,'position',sub_pos);
            hold on
            if c == 1 
                rc_num = 1;
                title_str = ['{\it',sprintf('%s}',freq_labels{f})];
            else
            end
            temp_data = [];
            temp_idx = [];
            sig_patch = [];
            for g = 1:length(subject_group)
                hold on
                marker_style = {'-o','LineWidth',l_width,'Color',group_colors(g,:),'markerfacecolor',[1 1 1],'MarkerSize',5};
                % compute error bars
                std_vals = nanstd(projected(g,q).rca(f).stats.SubjectAmpProjected(:,rc_num,:,c),0,3);
                err_vals = std_vals./sqrt(sum(~isnan(projected(g,q).rca(f).stats.SubjectAmpProjected(:,rc_num,:,c)),3));
                ci_vals = std_vals*1.96;
                [sig_vals, p_vals] = ttest(squeeze(projected(g,q).rca(f).stats.SubjectAmpProjected(:,rc_num,:,c)),0,'dim',2,'tail','right');
                
                for z = 1:length(sig_vals)
                    if sig_vals(z)
                        if z <= length(x_vals)
                            sig_patch_x = [x_vals(z)-log_step*0.53, x_vals(z)-log_step*0.53, x_vals(z)+log_step*0.53, x_vals(z)+log_step*0.53];
                            sig_patch_x(sig_patch_x < x_min) = x_min;
                            sig_patch_x(sig_patch_x > x_max) = x_max;
                            y_loc = y_max-(g-1)*y_unit/8;
                            sig_patch_y = [y_loc, y_loc-y_unit/8, y_loc-y_unit/8, y_loc];
                            sig_patch = [sig_patch, fill(sig_patch_x,sig_patch_y,group_colors(g,:),'edgecolor','none')];
                        else
                        end
                    else
                    end
                end
                
                
                % make line graphs
                if plot_medians
                    med_vals = median(projected(g,q).rca(f).stats.SubjectAmpProjected(:,rc_num,:,c),3);
                    h_val(g,1) = plot(x_vals,med_vals(1:10),marker_style{:});
                    if plot_binmean
                        h_val(g,2) = plot(extra_bins(g),med_vals(11),marker_style{:});
                        h_e = ErrorBars(...
                            [x_vals;extra_bins(g)],med_vals,[std_vals,std_vals], ...
                            'color',group_colors(g,:),'type','bar','cap',false,'barwidth',l_width);
                    else
                        h_e = ErrorBars(...
                            x_vals,med_vals,[std_vals(1:10),std_vals(1:10)], ...
                            'color',group_colors(g,:),'type','bar','cap',false,'barwidth',l_width);
                    end
                else
                    h_val(g,1) = plot(x_vals,projected(g,q).rca(f).stats.Amp(1:10,rc_num,c),marker_style{:});
                    if plot_binmean
                        h_val(g,2) = plot(extra_bins(g),projected(g,q).rca(f).stats.Amp(11,rc_num,c),marker_style{:});
                        h_e = ErrorBars(...
                            [x_vals;extra_bins(g)],projected(g,q).rca(f).stats.Amp(:,rc_num,c),[err_vals,err_vals], ...
                            'color',group_colors(g,:),'type','bar','cap',false,'barwidth',l_width);
                    else
                        h_e = ErrorBars(...
                            x_vals,projected(g,q).rca(f).stats.Amp(1:10,rc_num,c),[err_vals(1:10),err_vals(1:10)], ...
                            'color',group_colors(g,:),'type','bar','cap',false,'barwidth',l_width);
                    end   
                end
                arrayfun(@(x) uistack(x,'bottom'), h_val(g,:));
                cellfun(@(x) uistack(x,'bottom'), h_e);
                
                
                % get anova data
                temp_data = cat(1,temp_data, ...
                    permute(projected(g,q).rca(f).stats.SubjectAmpProjected(:,rc_num,:,c),[3,1,2]) );
                temp_idx = cat(1, temp_idx, ones(size(projected(g,q).rca(f).stats.SubjectAmpProjected(:,rc_num,:,c),3),1)*g );
            end
            hold on
            ref_y = plot([x_min, x_max], ones(2,1)*y_min, 'k-','linewidth',l_width);
            ref_x = plot(ones(2,1)*x_min, [y_min, y_max], 'k-','linewidth',l_width);
            uistack(ref_x,'bottom')
            uistack(ref_y,'bottom')
            arrayfun(@(x) uistack(x,'bottom'), sig_patch);
            
            anova_data(:,:,c+(q-1)*4) = temp_data;
            anova_groups(:,c+(q-1)*4) = temp_idx';
            hold on
            
            x_split = x_vals(end) + log_step*.5;
            if plot_binmean
                % plot dividing lines
                plot(ones(2,1)*x_split,[y_min,y_max],'k','LineWidth',l_width);
                plot(ones(2,1)*x_min,[y_min,y_max],'k','LineWidth',l_width);
                plot([x_min,x_max],ones(2,1)*y_min,'k','LineWidth',l_width);
            else
            end
            % plot noise patch, 
            % add mean bin noise twice to the end with extra bin
            % we are plotting the mean values over two x-axis points
            
            noise_set = max(cell2mat(arrayfun(@(x) projected(x,q).rca(f).stats.NoiseAmp(:,rc_num,c),1:length(subject_group),'uni',false)),[],2);
            %noise_set = mean(cell2mat(arrayfun(@(x) projected(x,q).rca(f).stats.NoiseAmp(:,rc_num,c),1:length(subject_group),'uni',false)),2);
            x_patch = [x_min,x_min,x_vals',x_split,x_split];
            y_patch = [0,noise_set(1),noise_set(1:end-1)',noise_set(end-1),0]; % start and end points just repeats of first and last
            %h_p(1) = fill(x_patch,y_patch,[.75 .75 .75],'edgecolor','none');
            if plot_binmean
                % additional vals for averages
                x_patch = [x_split,x_split,extra_bins,x_max,x_max];
                y_patch = [0,repmat(noise_set(end),1,5),0];
                %h_p(2) = fill(x_patch,y_patch,[.75 .75 .75],'edgecolor','none');
            else
            end
            %arrayfun(@(x) uistack(x,'bottom'), h_p);
            
            % plot layouts
            ylim([y_min,y_max]);
            xlim([x_min,x_max]);
            set(gca, gcaOpts{:},'ytick',y_min:y_unit:y_max,'xtick',xticks,'xticklabel',xtick_labels);
            if c == 1
                ylabel([title_str,': amplitude (\muV)'], text_params{:})
            elseif c == 4
                if ~mod(f,2)
                    [h_l,~,l_plots] = legend(h_val(:,1),{'typ','asd','adhd'}, text_params{:},'location','northeast','box','off');
                    for idx = 1:length(h_l.String)
                        h_l.String{idx} = ['\color[rgb]{' num2str(l_plots(idx).Color) '} ' h_l.String{idx}];
                    end
                    l_pos = get(h_l,'position');
                    l_pos(1) = l_pos(1) + l_pos(1) * 0.2;
                    l_pos(2) = l_pos(2) - l_pos(2) * 0.7;
                    set(h_l,'position',l_pos);
                end
            else
            end
            if mod(f,2)
                title(cond_names{c}, text_params{:});
            else
                if mod(c,2)
                    xlabel(x_str, text_params{:});
                else
                    xlabel(x_str, text_params{:});
                end
            end
            hold off
           
            
            % DO LME
            anova_ready = arrayfun(@(x) x, anova_data(:,1:10,c+(q-1)*4),'uni',false);
            group = repmat(anova_groups(:,c+(q-1)*4),length(bin_labels),1); 
            data = cat(2,anova_ready{:})';
            subject = repmat((1:length(ready_names))',length(bin_labels),1);
            bin = repmat(1:length(bin_labels),length(ready_names),1);
            bin = cellfun(@(x) str2double(x), bin_labels(bin(:)));
            group_lbl = arrayfun(@(x) sub_groups{x}, group,'uni',false);
            tbl = table(group_lbl,subject,bin,data);
            tbl.subject = categorical(tbl.subject);
            tbl.group_lbl = categorical(tbl.group_lbl,sub_groups);
            
            if fig_num == 1 || fig_num == 4
                if q == 1
                    r_table = sprintf('%s/analysis/r_vis_rc%d_%s_%d.csv',top_path,rc_num,freq_labels{f},c);
                else
                    r_table = sprintf('%s/analysis/r_aud_rc%d_%s_%d.csv',top_path,rc_num,freq_labels{f},c);
                end
                writetable(tbl, r_table, 'Delimiter',',','QuoteStrings',true)
            else
            end
            
            
            %model_out{c+(q-1)*4,f} = fitlme(tbl,'data ~ 1 + group_lbl*bin + (1|subject)','DummyVarCoding','reference');
            %anova_out{c+(q-1)*4,f} = anova(model_out{c+(q-1)*4,f},'DFMethod','satterthwaite');
            %if q == 1
            %    anova_path = sprintf('%s/analysis/anova_vis_rc%d_%s_c%d.mat',top_path,rc_num,freq_labels{f},c);
            %else
            %    anova_path = sprintf('%s/analysis/anova_aud_rc%d_%s_c%d.mat',top_path,rc_num,freq_labels{f},c);
            %end
            %save(anova_path,'data','subject','group','bin');
            
            % DO RANKSUM
            if any(isnan(anova_data))
                error('nans in anova_data');
            else
            end
            dset{1} = anova_data( anova_groups(:,c+(q-1)*4) == 1, :, c+(q-1)*4);
            dset{2} = anova_data( anova_groups(:,c+(q-1)*4) == 2, :, c+(q-1)*4);
            dset{3} = anova_data( anova_groups(:,c+(q-1)*4) == 3, :, c+(q-1)*4);
            
            % compute number of non-nan samples
            n = min(cell2mat(arrayfun(@(x) [sum(~isnan(dset{1}(:,x)));sum(~isnan(dset{2}(:,x)));sum(~isnan(dset{3}(:,x))) ], 1:length(bin_labels), 'uni',false)),[],2);
            n_str = arrayfun(@(x) sprintf('%d',x),n','uni',false);
            
            med_labels = cellfun(@(x,y) sprintf('%s m (n = %s)',x,y),sub_groups,n_str, 'uni',false );
            
            % compute medians
            str_median1 = [ med_labels{1}, arrayfun(@(x) num2str(nanmedian(dset{1}(:,x)),'%0.4f'),  1:(length(bin_labels)+1), 'uni',false) ];
            str_median2 = [ med_labels{2}, arrayfun(@(x) num2str(nanmedian(dset{2}(:,x)),'%0.4f'),  1:(length(bin_labels)+1), 'uni',false) ];
            str_median3 = [ med_labels{3}, arrayfun(@(x) num2str(nanmedian(dset{3}(:,x)),'%0.4f'),  1:(length(bin_labels)+1), 'uni',false) ];
             
            
            if q == 1
                table_path = sprintf('%s/analysis/vis_rc%d_%s.csv',top_path,rc_num,freq_labels{f});
                str_array = [str_array; sprintf('%s: %s',freq_labels{f},cond_names{c}),arrayfun(@(x) num2str(x,'bin%d'),1:length(bin_labels),'uni',false),'ave' ];
            else
                table_path = sprintf('%s/analysis/aud_rc%d_%s.csv',top_path,rc_num,freq_labels{f});
                str_array = [str_array; sprintf('%s: %s',freq_labels{f},cond_names{c}),arrayfun(@(x) num2str(x,'bin%d'),1:length(bin_labels),'uni',false),'ave' ];
            end
            
            str_array = [str_array; 'disp (arcmins)',bin_labels','-'];
            str_array = [str_array; str_median1; str_median2; str_median3]; 
            
            
            
            test_sets = [1,2; 1,3; 2,3];
            for t = 1:length(test_sets)
                label1 = sub_groups{test_sets(t,1)};
                label2 = sub_groups{test_sets(t,2)};
                % compute wilcoxon and p-values
                [w_p,w_h,w_stats] = arrayfun(@(x) ranksum(dset{test_sets(t,1)}(:,x),dset{test_sets(t,2)}(:,x),'method','approximate'),1:(length(bin_labels)+1));
                % compute t-tests
                % [t_h,t_p,~,t_stats] = arrayfun(@(x) ttest2(dset{test_sets(t,1)}(:,x),dset{test_sets(t,2)}(:,x),'vartype','unequal'),1:(length(bin_labels)+1),'uni',false)
                str_w = arrayfun(@(x) num2str(w_stats(x).ranksum,'%0.0f'),1:(length(bin_labels)+1),'uni',false);
                str_w = [sprintf('%s vs %s W',label1,label2),str_w];
                str_p = arrayfun(@(x) num2str(x,'%0.4f'),w_p,'uni',false);
                sig_idx = cell2mat(arrayfun(@(x) x < 0.0001,w_p,'uni',false));
                str_p(sig_idx) = {'<0.0001'};
                str_p = [sprintf('%s vs %s p',label1,label2),str_p];
                % compute AUROC
                mes_struct = arrayfun(@(x) mes(dset{test_sets(t,1)}(:,x),dset{test_sets(t,2)}(:,x),'auroc'),1:(length(bin_labels)+1));
                str_auroc = arrayfun(@(x) num2str(mes_struct(x).auroc,'%0.4f'),1:(length(bin_labels)+1),'uni',false);
                str_auroc = [sprintf('%s vs %s AUROC',label1,label2),str_auroc];
                str_array = [str_array; str_w; str_p; str_auroc];
            end
        end
        egi_h(fig_num) = subplot(2,5,[5,10]);
        hold on
        mrC.plotOnEgi(projected(1,q).rca(f).A(:,1));
        add_val = 2;
        new_pos = get(egi_h(fig_num) ,'position');
        new_pos(1) = new_pos(1)-(new_pos(3)*add_val/2);
        new_pos(2) = new_pos(2)-(new_pos(4)*add_val/2);
        new_pos(3:4) = new_pos(3:4)*(1+add_val);
        set(egi_h(fig_num),'position',new_pos);
        rc_pos = [0,1.75];%[max(get(gca,'xlim')),mean(get(gca,'ylim'))];
        if q == 1
            text(rc_pos(1),rc_pos(2), ['{\it', sprintf('visual RC no. %d}', rc_num)], text_params{:}, 'horizontalalignment', 'center');
        else
            text(rc_pos(1),rc_pos(2), ['{\it', sprintf('auditory RC no. %d}', rc_num)], text_params{:}, 'horizontalalignment', 'center');
        end
        % write table
        ready_table = array2table(str_array);
        clear str_*
        writetable(ready_table,table_path,'WriteRowNames',false,'WriteVariableNames',false);
    end
end
for f = 1:4
    if f == 1
        data_figpath = sprintf('%s/analysis/data_vis',top_path);
    elseif f == 4
        data_figpath = sprintf('%s/analysis/data_aud',top_path);
    else
        continue
    end
    figure(f)
    drawnow;
    set( figure(f), 'units', 'centimeters');
    fig_pos = get(figure(f),'pos');
    fig_pos(4) = 12;
    fig_pos(3) = 30;
    set(figure(f),'pos',fig_pos);
    export_fig(sprintf('%s.png', data_figpath),'-png','-opengl','-m5','-transparent',gcf);
end

%% COMPARISON PLOT
close all;
for q = 1:2
    figure;
    if q == 1
        comp_figpath = sprintf('%s/analysis/vis_comp',top_path);
        comp_conds = [1,3];
        comp_freqs = [1,2];
    else
        comp_figpath = sprintf('%s/analysis/aud_comp',top_path);
        comp_conds = [2,4];
        comp_freqs = [3,4];
    end
    for g = 1:3
        if g == 1
            if q == 1
                % x-values
                x_vals = cell2mat(cellfun(@(x) reallog(str2num(x)), vis_labels, 'uni',false));
                xtick_labels = [0.5,1,2,5,10,20,40,80];
                xticks = reallog(xtick_labels);
                bin_labels = vis_labels;
                x_str = 'contrast (%)';
            else
                x_vals = cell2mat(cellfun(@(x) reallog(str2num(x)), aud_labels, 'uni',false));
                x_vals = interp1([min(x_vals),max(x_vals)],[50,70],x_vals);
                xtick_labels = arrayfun(@(x) num2str(x,'%0.0f'), 50:5:70, 'uni',false);
                xticks = cellfun(@(x) str2num(x), xtick_labels);
                bin_labels = aud_labels;
                x_str = 'sound pressure (\it{dB})';
            end
            log_step = diff(x_vals(1:2)); % step size
            x_min = x_vals(1)-log_step*.5;
            if plot_binmean
                x_max = x_vals(end)+log_step*2.5; % add 2 steps
            else
                x_max = x_vals(end)+log_step*0.5;
            end
            extra_bins = arrayfun(@(x) x_vals(end)+x, [log_step,log_step*1.5,log_step*2]);
        else
        end
        hold on
        for f = 1:length(comp_freqs)
            if q == 1
                if f ==1
                    y_min = 0;
                    y_unit = 4;
                    y_max = 16;
                else
                    y_min = 0;
                    y_unit = 2;
                    y_max = 8;
                end
            else
                if f ==1
                    y_min = 0;
                    y_unit = 2;
                    y_max = 8;
                else
                    y_min = 0;
                    y_unit = 1;
                    y_max = 4;
                end
            end
            title_str = ['{\it',sprintf('%s}',freq_labels{comp_freqs(f)})];
            subplot(2,3,g+(f-1)*3);
            for c = 1:length(comp_conds)
                if f == 1
                    m_color = [1 1 1];
                else
                    m_color = [1,1,1]; %group_colors(g,:);
                end
                if c == 1
                    marker_style = {'-o','LineWidth',l_width,'Color',group_colors(g,:),'markerfacecolor', m_color ,'MarkerSize',6};
                else
                    marker_style = {'-d','LineWidth',l_width,'Color',group_colors(g,:),'markerfacecolor', m_color,'MarkerSize',6};
                end
                std_vals = nanstd(projected(g,q).rca(f).stats.SubjectAmpProjected(:,rc_num,:,comp_conds(c)),0,3);
                err_vals = std_vals./sqrt(sum(~isnan(projected(g,q).rca(comp_freqs(f)).stats.SubjectAmpProjected(:,rc_num,:,comp_conds(c))),3));
                ci_vals = std_vals*1.96;
                if plot_medians
                    med_vals = median(projected(g,q).rca(f).stats.SubjectAmpProjected(:,rc_num,:,comp_conds(c)),3);
                    h_val(g,1) = plot(x_vals,med_vals(1:10),marker_style{:});
                    hold on
                    if plot_binmean
                        h_val(g,2) = plot(extra_bins(g),med_vals(11),marker_style{:});
                        h_e = ErrorBars(...
                            [x_vals;extra_bins(g)],med_vals,[std_vals,std_vals], ...
                            'color',group_colors(g,:),'type','bar','cap',false,'barwidth',l_width);
                    else
                        h_e = ErrorBars(...
                            x_vals,med_vals,[std_vals(1:10),std_vals(1:10)], ...
                            'color',group_colors(g,:),'type','bar','cap',false,'barwidth',l_width);
                    end
                else
                    h_val(g,1) = plot(x_vals,projected(g,q).rca(comp_freqs(f)).stats.Amp(1:10,rc_num,comp_conds(c)),marker_style{:});
                    hold on
                    if plot_binmean
                        h_val(g,2) = plot(extra_bins(g),projected(g,q).rca(comp_freqs(f)).stats.Amp(11,rc_num,comp_conds(c)),marker_style{:});
                        h_e = ErrorBars(...
                            [x_vals;extra_bins(g)],projected(g,q).rca(comp_freqs(f)).stats.Amp(:,rc_num,comp_conds(c)),[err_vals,err_vals], ...
                            'color',group_colors(g,:),'type','bar','cap',false,'barwidth',l_width);
                    else
                        h_e = ErrorBars(...
                            x_vals,projected(g,q).rca(comp_freqs(f)).stats.Amp(1:10,rc_num,comp_conds(c)),[err_vals(1:10),err_vals(1:10)], ...
                            'color',group_colors(g,:),'type','bar','cap',false,'barwidth',l_width);
                    end   
                end
                arrayfun(@(x) uistack(x,'top'), h_val(g,:));
                cellfun(@(x) uistack(x,'bottom'), h_e);
                hold on
            end
            % plot layouts
            ylim([y_min,y_max]);
            xlim([x_min,x_max]);
            set(gca, gcaOpts{:},'ytick',y_min:y_unit:y_max,'xtick',xticks,'xticklabel',xtick_labels);
            if f == 1
                title(sub_groups{g}, text_params{:}); 
                if g == 1
                    leg_p(1) = plot(nan, nan, '-ko', 'LineWidth', l_width, 'markerfacecolor', [1,1,1] ,'MarkerSize', 6);
                    leg_p(2) = plot(nan, nan, '-kd', 'LineWidth', l_width, 'markerfacecolor', [1,1,1] ,'MarkerSize', 6);
                    legend(leg_p,{'unimodal','crossmodal'}, 'box','off',text_params{:}, 'location', 'northwest');
                else
                end
            else
                if g == 1
                    xlabel(x_str, text_params{:});
                else
                end
            end
            if g == 1
                ylabel([title_str,': amplitude (\muV)'], text_params{:})
            else
            end
        end
        hold off
    end
    set( gcf, 'units', 'centimeters');
    fig_pos = get(gcf,'pos');
    fig_pos(4) = 12;
    fig_pos(3) = 24;
    set(gcf,'pos',fig_pos);
    export_fig(sprintf('%s.png', comp_figpath),'-png','-opengl','-m5','-transparent',gcf);
end