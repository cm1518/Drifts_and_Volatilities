%% ========================================================================
% GIRF Figures in the paper
% -------------------------------------------------------------------------
%
% Figure 5: True vs. Observed data
% Figure 6: Evolving volatilities
% Figure 7: (Forecast-Based)Conditional Mean
%
%==========================================================================

clear('all'); close('all'); clc;

% Make sure to have the same cd as the main script <Main_me_tvpvar_public.m>
main_path = [pwd '\'];
save_path = [pwd '\Results\'];

addpath(main_path)
addpath(save_path)
addpath([pwd '\Data\'])
addpath([pwd '\MiscCodes\'])

cd(main_path)

prior_sel   = 1;
data_sel    = 2;

switch data_sel
    case 1
        var_name    = {'Real GDP','Inflation','Interest Rate','Spread','M2 Growth'};
    case 2
        var_name    = {'Real GDP','Inflation','Interest Rate','Spread','Money Growth'};
end

varname     = var_name;
save_name   = 'FinalRevisionQE';


%% I. Load necessary Mat files
% =========================================================================
str_small_girf =  '2AAMW_ME_Level_Model_2_Prior_1_Yres_1_SRHorz_3_sr_GIRF';
str_large_girf =  'AAMW_ME_Level_Model_2_Prior_1_Yres_1_SRHorz_3_sr_GIRF';

eval(sprintf('load %s',str_large_girf))
% eval(sprintf('load AAMW_ME_Level_%d_Prior_%d_Yres_%d_SRHorz_%d_sr_GIRF.mat',data_sel,prior_sel,DO_Y,sr_horz))

% -------------------------------------------------------------------------
% Figure settings
% -------------------------------------------------------------------------
year_indx   = [find(yearlab==1931),find(yearlab==1937),find(yearlab==1945),...
    find(yearlab==1965),find(yearlab==1974),find(yearlab==1987),...
    find(yearlab==1996),find(yearlab==2001),find(yearlab==2008)];

RGB         = [1 1 1];
DO_SAVE     = 1;


%% Historical Charimanship Figures
% =========================================================================
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

varname  = {sprintf('GDP Growth'),sprintf('Inflation'),sprintf('Interest Rate'),...
    sprintf('Spread'),sprintf('Money Growth')};

color_order = [ 0    0  .54;
    .54   0   0;
    0    0   0;
    0  .54 .54;
    .54   0 .54];

rgb_order   = [[173	216	230]/255;
    [255	160	122]/255;
    .7  .7  .7;
    [193 255 193]/255;
    [255	187	255]/255];

%% Figure 5 and Figure 7
% =========================================================================
for vers_i = 1 : 2
    switch vers_i
        case 1
            tmp_irf = tube.srgirf; save_name = 'HistoricalSRGIRF';
        case 2
            tmp_irf = tube.norm; save_name = 'HistoricalSRGIRF_NormInterestRate';
    end
    
    tmpTime = [1:nhorz nhorz:-1:1]';
    sp_indx = 0;
    clear hh
    FontSize = 4;
    for ii = 1 : M
        ymin1 = []; ymax1 = [];
        for jj = 1:7;
            clear alpha
            sp_indx = sp_indx + 1;
            hh(sp_indx)=subplot(5,7,sp_indx);hold('on');
            eval(sprintf('time_tmp=time_index%s;',int2str(jj)))
            tmp_aa = squeeze(tmp_irf(ii,time_tmp,1:end,:));
            tmp_aa = permute(tmp_aa,[2 1 3]);
            tmp_aa = tmp_aa(:,:);
            tmp_aa = squeeze(prctile(tmp_aa,[16 50 84],2));
            
            if ~any(ii==[4]);
                plot(zeros(1,nhorz)+sr_horz,100*(-floor(nhorz/2):floor(nhorz/2)),'k-','LineWidth',.5);
            end
            plot(1:nhorz,zeros(1,nhorz),'k-','LineWidth',.5);
            tmpGIRF = [tmp_aa(:,1);flipdim(tmp_aa(:,3),1)];
            patch(tmpTime,tmpGIRF,rgb_order(ii,:),'EdgeColor',color_order(ii,:))
            plot(1:nhorz,tmp_aa(:,2),'Color',color_order(ii,:),'LineWidth',1.5);
            axis tight,grid on;
            
            ymin1 = [ymin1;floor(min(tmp_aa(:)*10))/10];
            ymax1 = [ymax1; ceil(max(tmp_aa(:)*10))/10];
            
            
            if sp_indx < 8;title(tit_name(sp_indx),'FontSize',FontSize,'FontWeight','Bold','FontName','Calibri');end
            if sp_indx == 1;     ylabel(varname(1),'FontSize',FontSize,'Color','k','FontWeight','Bold','FontName','Calibri');
            elseif sp_indx == 8; ylabel(varname(2),'FontSize',FontSize,'Color','k','FontWeight','Bold','FontName','Calibri');
            elseif sp_indx == 15;ylabel(varname(3),'FontSize',FontSize,'Color','k','FontWeight','Bold','FontName','Calibri');
            elseif sp_indx == 22;ylabel(varname(4),'FontSize',FontSize,'Color','k','FontWeight','Bold','FontName','Calibri');
            elseif sp_indx == 29;ylabel(varname(5),'FontSize',FontSize,'Color','k','FontWeight','Bold','FontName','Calibri');
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
            set(gca,'XColor',[.5 .5 .5]-.3,'YColor',[.5 .5 .5]-.3,'box','off','LineWidth',.5,'GridLineStyle',':')
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
    print(gcf,'-dtiff','-r1200',sprintf('%s%s_Model_%d_Prior_%d',save_path,save_name,data_sel,prior_sel))
    close all
end

%% Figure 6
% note that for figure 6 takes some time!
do_figure_6_public



