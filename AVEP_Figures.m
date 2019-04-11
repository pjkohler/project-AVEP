%% PLOT RCA TRAINED ON EACH CONDITION SEPARATELY
close all;

codeFolder = '/Users/kohler/code';
rcaCodePath = sprintf('%s/git/rcaBase', codeFolder);
addpath(genpath(rcaCodePath));
addpath(genpath(sprintf('%s/git/mrC', codeFolder)));
addpath(genpath(sprintf('%s/git/schlegel/matlab_lib', codeFolder)));
setenv('DYLD_LIBRARY_PATH','')

% codeFolder = '/Users/labmanager/Desktop/LabManager/MatAnal';
% rcaCodePath = sprintf('%s/rcaBase',codeFolder);
% addpath(genpath(rcaCodePath));
% addpath(genpath(sprintf('%s/mrC',codeFolder)));
% addpath(genpath(sprintf('%s/mrC/+mrC',codeFolder)));
% addpath(genpath(sprintf('%s/schlegel/matlab_lib',codeFolder)));
% addpath(genpath(sprintf('%s/matlab_lib',codeFolder)));
% addpath(genpath(sprintf('%s/export_fig',codeFolder)));
% setenv('DYLD_LIBRARY_PATH','')



%load('AVEPOutputf1F_4_20180501')
load('AVEPOutputf5F_8_20180508')

for c=1:(length(condNames)-1)    
    % set figure size in the beginning
    figure;
    set(gcf, 'units', 'centimeters');
    figPos = get(gcf,'pos');
    figPos(3) = 5*(length(freqsToUse)+1);
    figPos(4) = 30;
    set(gcf,'pos',figPos);
    plotSettings.titleOn = true;
    plotSettings.titleToUse = [];
    legendLabels = {sprintf('ASD,%s  %s  n= %d', suffix, suffix, length(ASDpathnames)),sprintf('TYP,%s  %s  n= %d', suffix, suffix, length(TYPpathnames)),...
        sprintf('ADHD, n= %d', length(ADHDpathnames))};
    
    if strfind(condNames{c}, 'VIS')
        plotSettings.xlabel = 'Contrast';
        plotSettings.xTick = [0, .1, 1, 10, 100];
    elseif strfind(condNames{c},'AUD')
        plotSettings.xlabel = 'Loudness';
        plotSettings.xTick = [0, .01, .05, .1, 5];
    end
    
    %plotSettings.xTick = fullASD.(condNames{c}).data.settings.binLevels{1}([1 5 10]);
    %plotSettings.xTick = [0, .1, 1, 10, 100];
    binVals = fullASD.(condNames{c}).settings.binLevels{1};
    if iscell(binVals)
    binVals = binVals';    %work around for binVals being a cell array with str variables. transpose, convert to double, convert binVals to vector.
    for binIndex = 1:length(binVals)
        binVals{binIndex} = str2double(binVals{binIndex});
    end
    binVals = cell2mat(binVals);
    end
    
    nRC = fullASD.(condNames{c}).settings.nComp;
    nFreq = length(freqsToUse);
    plotSettings.ymax = [];

    %cell 2 vector
    %rcaSettings.binLevels = [rcaSettings.binLevels{:}];
    
    % plot settings
    lWidth = 1;
    errorWidth = 0.3;
    fSize = 10;
    cBrewer = load('colorBrewer.mat');
    subColors = [cBrewer.rgb20(1,:); cBrewer.rgb20(5,:); cBrewer.rgb10(2,:)];
    gcaOpts = {'tickdir','out','ticklength',[0.0500,0.0500],'XTick', plotSettings.xTick, 'XScale', 'log', 'box','off','fontsize',fSize,'fontname','Arial','linewidth',lWidth};
    yMax = [9,3,3,3,3,9;
             4,2,2,2,2,4;
             4,2,2,2,2,4;
             2,.5,.5,.5,.5,2;
             .4,.2,.2,.2,.2,.4;
             .4,.2,.2,.2,.2,.4];

    for r = 1:(nRC+1)  
        if r<6
            egiH(r) = subplot(nRC+1,nFreq+1,(r-1)*nFreq+nFreq+r);
            hold on
            rcaColorBar = [min(fullASD.(condNames{c}).A(:,r)),max(fullASD.(condNames{c}).A(:,r))];
            newExtreme = max(abs(rcaColorBar));
            rcaColorBar = [-newExtreme,newExtreme];
            mrC.plotOnEgi(fullASD.(condNames{c}).A(:,r),rcaColorBar);
            hold off
        else
        end

        for f=1:nFreq
            subplot(nRC+1,nFreq+1,(r-1)*nFreq+r+(f-1));
            hold on
            for s=2
                     ampH(s)= errorbar(binVals,ampVals.(condNames{c})(:,f,r,s),errLB.(condNames{c})(:,f,r,s),errUB.(condNames{c})(:,f,r,s),'o-','CapSize',0, 'MarkerSize', 2,'Color', subColors(s,:), 'MarkerEdgeColor', subColors(s,:), 'MarkerFaceColor', subColors(s,:),'LineWidth',errorWidth);
                     %ampH(s)= semilogx(binVals,ampVals.(condNames{c})(:,f,r,s),'o-', 'MarkerSize', 2,'Color', subColors(s,:), 'MarkerEdgeColor', subColors(s,:), 'MarkerFaceColor', subColors(s,:),'LineWidth',errorWidth);
                
%                  if ~isnan(NR_pOpt.(condNames{c})(1,f,r,s))
%                      % plot Naka-Rushton
%                      nFine = 1e2;
%                      nrX = linspace( min(binVals), max(binVals), nFine )';
%                      nrVals = hModel( nrX, NR_pOpt.(condNames{c})(:,f,r,s));
%                      plot( nrX, nrVals, '-','Color', subColors(s,:), 'LineWidth',lWidth);
%                  else
%                  end
            end
            set(gca,gcaOpts{:}, 'XScale','log');%,'ytick',linspace(0,yMax(f,r),3)); 
            xlim([(min(binVals)*.75), (max(binVals)+max(binVals)*0.3)]);
            %ylim([0,yMax(f,r)]);
            if f==1
                ylabel('Amplitude (\muVolts)')
                if r > nRC
                    text(min(binVals)*.01,max(get(gca,'ylim'))*1.05,'Ch. 4','fontsize',fSize,'fontname','Arial');
                else
                    text(min(binVals)*.01,max(get(gca,'ylim'))*1.05,['RC ',num2str(r)],'fontsize',fSize,'fontname','Arial');
                end
            else
            end
            if r==1
                title(sprintf('%s: %s',condNames{c}, fullASD.(condNames{c}).settings.freqLabels{1}{f}),'fontsize',fSize,'fontname','Arial');
            elseif r==6 
                xlabel(plotSettings.xlabel);
                if f==nFreq
                   % hLeg = legend(ampH,{'ASD','TYP','ADHD'},'location','northeast');
                     hLeg = legend(ampH,legendLabels,'location','northeast');
                      lPos = get(hLeg,'pos');
                      lPos(1) = lPos(1)+lPos(1)*.3;   %was .2
                      set(hLeg,'pos',lPos);
                else
                end
            end
            hold off;
        end
    end
    drawnow;
    for r = 1:6
        if r < 6
            addVal = 1/2;
            newPos = get(egiH(r),'position');
            newPos(1) = newPos(1)-(newPos(3)*addVal/2);
            newPos(2) = newPos(2)-(newPos(4)*addVal/2);
            newPos(3:4) = newPos(3:4)*(1+addVal);
            set(egiH(r),'position',newPos);
        else
             lPos = get(hLeg,'pos');
             lPos(1) = lPos(1)+lPos(1)*.03;
             set(hLeg,'pos',lPos);
        end
            
    end

        export_fig (sprintf('201804_AVEP_%s_%dF%dF',condNames{c}, min(freqsToUse), max(freqsToUse)), '-pdf', '-eps', '-transparent', gcf);
 
end
