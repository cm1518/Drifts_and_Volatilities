eval(sprintf('load %sAAMW_ME_Level_Model_%d_Prior_%d',save_path,data_sel,prior_sel))

number_of_draws = 1;        % You evaluate on st mean or meadia, hence 1!
nalpha          = 1000;
NMCInteg        = nalpha;
SC              = 0;
final_kept_draw = number_of_draws*nalpha;

nhorz           = 21;
h               = nhorz;
sr_horz         = 3;
DO_Y            = 1;
sr_end          = zeros(number_of_draws,t);
draw_index      = randsample(nrep,number_of_draws);

FullSet.srgirf     = nan(M,t,nhorz,final_kept_draw);

save_index      = [1:nalpha:final_kept_draw final_kept_draw+1];
save_index1     = [1:3:number_of_draws*3 number_of_draws*3+1];

% Choose to either evaluate get the full identfied set at the posteior
% meanor medain
% -------------------------------------------------------------------------
% Bt_med      = median(Bt_post,3);
% At_med      = median(At_post,3);
% St_med      = median(Sigt_post,3);
% sig_med     = median(sig_post,3);

Bt_mean     = mean(Bt_post,3);
At_mean     = mean(At_post,3);
St_mean     = mean(Sigt_post,3);
sig_mean    = mean(sig_post,3);

tic;

for l = 1 : number_of_draws
    for i = 1 : t
        if mod(i,1)==0;
            fprintf('Rotation for Time\t%3.0f\t(%4.2f\t minutes.)\n',i,toc/60);
        end
        g_rand = l;
        %% Here we do GIRFs á la Koop, Pesaran and Potter (1996)
        % =====================================================================
        Btmp    = Bt_mean(:,i);
        Atmp    = At_mean(:,i);
        Sigtmp  = St_mean(:,i);
        Qtmp    = squeeze(mean(Q_post,3));
        Wtmp    = squeeze(mean(W_post,3));
        Stmp    = squeeze(mean(S_post,3));
        stem    = diag(sig_mean(i,:));
        
        nhor    = 21;
        tt      = p+1;
        YYinit  = Y(length(Y)-t+1-p:length(Y),:)';
        if i <= p
            YY      = YYinit(:,i:i+p-1,:);
        else
            YY      = y(:,i-p:i-1);
        end
        [ind,girf,sirf]     = do_sr_girf_gr(NMCInteg,SC,K,M,numa,nhor,p,Btmp,Atmp,Sigtmp,Qtmp,Wtmp,Stmp,YY,DO_Y,sr_horz);
        
        FullSet.srgirf(:,i,:,save_index(1):save_index(2)-1) = permute(sirf,[2 3 1]);
    end
    save_index(1)=[];
    save_index1(1)=[];
end
fprintf('\n')
fprintf('==================================================================\n')
fprintf('Full GIRF Run: %8.0f\t(%4.2f\t minutes.)\n',l,toc/60)
fprintf('==================================================================\n')

%% Do Figures of identfied set at median
year_indx   = [find(yearlab==1931),find(yearlab==1937),find(yearlab==1945),...
    find(yearlab==1965),find(yearlab==1974),find(yearlab==1987),...
    find(yearlab==1996),find(yearlab==2001),find(yearlab==2008)];

RGB         = [1 1 1];
DO_SAVE     = 1;
nhor =nhorz;
time_index1 = 1:find(yearlab==1929);                        % Fed's Formative Years
time_index2 = length(time_index1)+1:find(yearlab==1941);    % Great Depression
time_index3 = length(time_index2)+1:find(yearlab==1951);    % WW-II and aftermath
time_index4 = length(time_index3)+1:find(yearlab==1965);    % Treasury Fed Accord to mid-60
time_index5 = length(time_index4)+1:find(yearlab==1982);    % Great Inflation
time_index6 = length(time_index5)+1:find(yearlab==2007);    % Great Moderation
time_index7 = length(time_index6)+1:t;                      % Great Recession and Aftermath

tit_name = {sprintf('\nFed''s Formative Years\n\n\t1913 - 1929'),...
    sprintf('\nGreat Depression\n\n\t1929 - 1941'),...
    sprintf('\nWorld War II \nand  Aftermath\n\t1941 - 1951'),...
    sprintf('\nTreasury Fed Accord\n\n\t1951 - 1965'),...
    sprintf('\nGreat Inflation\n\n\t1965 - 1982'),...
    sprintf('\nGreat Moderation\n\n\t1982 - 2007'),...
    sprintf('\nGreat Recession\nand  Aftermath\n\t2007 - 2012')};

color_order = [ 0    0  .54;
    .54   0   0;
    0    0   0;
    0  .54 .54;
    .54   0 .54];

rgb_order   = [[173	216	230]/255;
    [255	160	122]/255;
    .7  .7  .7;
    [155,205,155]/255;
    [255	187	255]/255];

rgb_order2  = [
    [224,255,255]/255;
    [255,218,185]/255;
    .9  .9  .9;
    [193 255 193]/255;
    [255,225,255]/255];

%% Figure 6
% =========================================================================
eval(sprintf('load AAMW_ME_Level_Model_%d_Prior_%d_Yres_%d_SRHorz_%d_sr_GIRF.mat',data_sel,prior_sel,DO_Y,sr_horz))

save_name   = 'HistoricalPosterior90vsIdentifiedSet';
tmp_irf     = tube.srgirf; 
tmp_firf    = FullSet.srgirf;

tmpTime = [1:nhor nhor:-1:1]';
sp_indx = 0;
clear hh
FontSize = 4;
for ii = 1 : M
    ymin1 = []; ymax1 = [];
    for jj = 1:7;
        clear alpha
        sp_indx = sp_indx + 1;
        hh(sp_indx)=subplot(5,7,sp_indx);
        hold('on');
        
        eval(sprintf('time_tmp=time_index%s;',int2str(jj)))
        
        % Get posterior 90% 
        tmp_aa = squeeze(tmp_irf(ii,time_tmp,1:end,:));
        tmp_aa = permute(tmp_aa,[2 1 3]);
        tmp_aa = tmp_aa(:,:);
        tmp_aa = squeeze(prctile(tmp_aa,[5 95],2));
        
        % Get full identified set at posterior median/mean
        tmp_faa = squeeze(tmp_firf(ii,time_tmp,1:end,:));
        tmp_faa = permute(tmp_faa,[2 1 3]);
        tmp_faa = tmp_faa(:,:);
        tmp_faa = squeeze(prctile(tmp_faa,[0 100],2));
        
        if ~any(ii==[4]);
            plot(zeros(1,nhor)+sr_horz,100*(-floor(nhor/2):floor(nhor/2)),'k-','LineWidth',.5);
        end
        plot(1:nhor,zeros(1,nhor),'k-','LineWidth',.5);
        
        tmpGIRF         = [tmp_aa(:,1);flipdim(tmp_aa(:,2),1)];
        tmpGIRF_FullSet = [tmp_faa(:,1);flipdim(tmp_faa(:,2),1)];
        
        % Posterior percentile (90%)
        % -------------------------------------------------------------
        patch(tmpTime,tmpGIRF,color_order(ii,:),'EdgeColor',color_order(ii,:))
        
        % Identified set at posterior Mean
        % -------------------------------------------------------------
        patch(tmpTime,tmpGIRF_FullSet,rgb_order2(ii,:),'EdgeColor',rgb_order(ii,:))
        
        axis tight,grid on;
        
        ymin1 = [ymin1;floor(min([tmp_aa(:);tmp_faa(:)]*10))/10];
        ymax1 = [ymax1; ceil(max([tmp_aa(:);tmp_faa(:)]*10))/10];
        
        if sp_indx < 8;title(tit_name(sp_indx),'FontSize',FontSize,'FontWeight','Bold','FontName','Calibri');end
        if sp_indx == 1;     ylabel(var_name(1),'FontSize',FontSize,'Color','k','FontWeight','Bold','FontName','Calibri');
        elseif sp_indx == 8; ylabel(var_name(2),'FontSize',FontSize,'Color','k','FontWeight','Bold','FontName','Calibri');
        elseif sp_indx == 15;ylabel(var_name(3),'FontSize',FontSize,'Color','k','FontWeight','Bold','FontName','Calibri');
        elseif sp_indx == 22;ylabel(var_name(4),'FontSize',FontSize,'Color','k','FontWeight','Bold','FontName','Calibri');
        elseif sp_indx == 29;ylabel(var_name(5),'FontSize',FontSize,'Color','k','FontWeight','Bold','FontName','Calibri');
        end
        if sp_indx == 32; xlabel(sprintf('Horizon in Years\n'),'FontSize',FontSize,'FontWeight','Bold','Color','k','FontName','Calibri'); end
        set(gca,'XTick',[5:4:21]);
        if sp_indx <=28;set(gca,'XTickLabel',[]);
        else
            set(gca,'XTick',[5:4:21],'XTickLabel',{'1','2','3','4','5'},'FontWeight','Bold','FontName','Calibri');
        end
        if ~any(sp_indx == [1:7:M*7])
            set(gca,'YTickLabel',[],'FontWeight','Bold');
        end
        
        set(gca,'FontSize',FontSize-2)
        set(gca,'XColor',[.5 .5 .5]-.3,'YColor',[.5 .5 .5]-.3,'box','off',...
            'LineWidth',.5,'GridLineStyle',':')%,'FontWeight','Bold')
        alpha(.5)
    end
    ymin(ii) = floor(min(ymin1(:)*10))/10;
    ymax(ii) =  ceil(max(ymax1(:)*10))/10;
end
set(hh( 1: 7),'YLim',[ymin(1) ymax(1)],'FontWeight','Bold')
set(hh( 8:14),'YLim',[ymin(2) ymax(2)],'FontWeight','Bold')
set(hh(15:21),'YLim',[ymin(3) ymax(3)],'FontWeight','Bold')
set(hh(22:28),'YLim',[ymin(4) ymax(4)],'FontWeight','Bold')
set(hh(29:35),'YLim',[ymin(5) ymax(5)],'FontWeight','Bold')
set(hh,'FontSize',FontSize,'FontWeight','Bold')

saveTightFigure(gcf,sprintf('%s%s_Model_%d_Prior_%d.tiff',save_path,save_name,data_sel,prior_sel))
delete(sprintf('%s%s_Model_%d_Prior_%d.tiff',save_path,save_name,data_sel,prior_sel))
print(gcf,'-dtiff','-r1200',sprintf('_Model_%d_Prior_%d',save_path,save_name,data_sel,prior_sel))

close all

