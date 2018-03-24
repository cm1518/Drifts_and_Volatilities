%% ========================================================================
% Descriptive Figures in the paper
% -------------------------------------------------------------------------
% 
% Figure 1: True vs. Observed data
% Figure 2: Evolving volatilities
% Figure 3: (Forecast-Based)Conditional Mean
% Figure 4: (Forecast-Based) Correlation 
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
%eval(sprintf('load %sAAMW_ME_Level_Model_%d_Prior_%d',save_path,data_sel,prior_sel));

% -------------------------------------------------------------------------
% Figure 1: True vs. Observed data
% -------------------------------------------------------------------------
yearlabn        = yearlab(5:end); 
me_post_diff    = zeros(size(me_post));
var_i           = 1;
beg_me          = [5,5,1,1,5];          % Output growth, inflation, intereste rate, spread, money growth
end_me          = [130,130,22,22,178];  % variable corresponding last Mreak date for ME 
                                        % elements to be selected from yearlab

for var_i = 1 : M
    if beg_me(var_i)==5
        me_post_diff(beg_me(var_i):end_me(var_i),var_i,:) = me_post(beg_me(var_i):end_me(var_i),var_i,:) - me_post(1:end_me(var_i)-(beg_me(var_i)-1),var_i,:);
    else
        me_post_diff(beg_me(var_i):end_me(var_i),var_i,:) = me_post(beg_me(var_i):end_me(var_i),var_i,:);% Start at mal_me_lag+1==5 to align begining date of gr ME
    end
end

tmp_me          = repmat(y_obs',[1 1 nrep])-me_post_diff;
tmp_me(1:4,:,:) = [];
y_obs_tmp       = [y_obs(:,5:end)];
tmp_ax          = [yearlabn;flipdim(yearlabn,1)];
tmp_mep         = prctile(tmp_me,[16 50 84],3);
RGB_tmp         = [32 178 170]/255;
FontSize        = 4;
for ii = 1 : M+1
    subplot(3,2,ii);
    if ii>M
        hold('on')
        ggg(1)=patch(tmp_ax,tmp_i,RGB_tmp,'EdgeColor',RGB_tmp);
        ggg(2)=plot(yearlabn,squeeze(tmp_mep(:,ii-1,2)),'Color',[0 128 128]/255);             % ME Median
        ggg(3)=plot(yearlab(end_me(ii-1)+1:end),y_obs(ii-1,end_me(ii-1)+1:end),'Color',[0 0 0]);
        ggg(4)=plot(yearlab(5:end_me(ii-1)),y_obs(ii-1,5:end_me(ii-1)),'Color',[1 0 0],'LineStyle','-');
        ggg(5)=line([yearlab(find(yearlab==1930));yearlab(find(yearlab==1930))],[min_tmp; max_tmp],'LineWidth',1,'Color',[47 79 79]/255);
        alpha(.5)
        legend(ggg,sprintf('"True" Data (68%% Posterior)'),sprintf('"True" Data (Median)'),'Observed Data (w/o ME)','Observed Data (w/ ME)','Break Points in Measurement')
        title('Legend','FontSize',FontSize,'FontWeight','Bold')
        set(gca,'FontSize',FontSize,'FontWeight','Bold','Xlim',[0 1],'Ylim',[0 1],'XTickLabel',[],'YTickLabel',[],'Box','On')
    else
        hold('on')
        tmp_i = [squeeze(tmp_mep(:,ii,1));flipdim(squeeze(tmp_mep(:,ii,3)),1)];
        patch(tmp_ax,tmp_i,RGB_tmp,'EdgeColor',RGB_tmp);
        plot(yearlabn,squeeze(tmp_mep(:,ii,2)),'Color',[0 128 128]/255);
        plot(yearlab(end_me(ii)+1:end),y_obs(ii,end_me(ii)+1:end),'Color',[0 0 0]);
        plot(yearlab(5:end_me(ii)),y_obs(ii,5:end_me(ii)),'Color',[1 0 0],'LineStyle','-');
        axis('tight');
        grid('on')
        set(gca,'FontSize',FontSize,'Xlim',[1914 1959])
        min_tmp = floor(min([y_obs(ii,:) min(min(tmp_mep(:,ii,:)))])*10)/10;
        max_tmp = ceil( max([y_obs(ii,:) max(max(tmp_mep(:,ii,:)))])*10)/10;
        alpha(.5)
        if ii==1;
            line([yearlab(find(yearlab==1930));yearlab(find(yearlab==1930))],[min_tmp; max_tmp],'LineWidth',1,'Color',[47 79 79]/255)
            line([yearlab(find(yearlab==1947)-1);yearlab(find(yearlab==1947)-1)],[min_tmp; max_tmp],'LineWidth',1,'Color',[47 79 79]/255)
        elseif ii==2;
            line([yearlab(find(yearlab==1947)-1);yearlab(find(yearlab==1947)-1)],[min_tmp; max_tmp],'LineWidth',1,'Color',[47 79 79]/255)
        elseif ii==3;
            line([yearlab(find(yearlab==1920)-1);yearlab(find(yearlab==1920)-1)],[min_tmp; max_tmp],'LineWidth',1,'Color',[47 79 79]/255)
        elseif ii==4;
            line([yearlab(find(yearlab==1920)-1);yearlab(find(yearlab==1920)-1)],[min_tmp; max_tmp],'LineWidth',1,'Color',[47 79 79]/255)
        elseif ii==5;
            line([yearlab(find(yearlab==1918));yearlab(find(yearlab==1918))],[min_tmp; max_tmp],'LineWidth',1,'Color',[47 79 79]/255)
            line([yearlab(find(yearlab==1936));yearlab(find(yearlab==1936))],[min_tmp; max_tmp],'LineWidth',1,'Color',[47 79 79]/255)
            line([yearlab(find(yearlab==1959)-1);yearlab(find(yearlab==1959)-1)],[min_tmp; max_tmp],'LineWidth',1,'Color',[47 79 79]/255)
        end
        set(gca,'FontSize',FontSize,'Xlim',[1914 1960],'YLim',[min_tmp max_tmp])
        title(sprintf('%s',char(var_name(ii))),'FontSize',FontSize,'FontWeight','Bold')
    end
end
saveTightFigure(gcf,sprintf('%sME_Final_%d_Prior_%d.tiff',save_path,data_sel,prior_sel))
delete(sprintf('%sME_Final_%d_Prior_%d.tiff',save_path,data_sel,prior_sel))
print(gcf,'-dtiff','-r1200',sprintf('%sME_Final_%d_Prior_%d',save_path,data_sel,prior_sel))
close all

% -------------------------------------------------------------------------
% Figure 2: Evolving volatilities
% -------------------------------------------------------------------------
Fontsize = 4;
clear ggg
tmpa = [yearlab;flipdim(yearlab,1)];
figure('Name',sprintf('Evolving Volatility: Model %d and Prior %d',data_sel,prior_sel));
for ii = 1 : size(sig_post,2)+1
    subplot(3,2,ii);hold('on')
    if ii>M
        hold('on')
        ggg(2)=patch(tmpa,tmpb,[255 182 193]./255,'EdgeColor',[255 182 193]./255);
        ggg(1)=plot(yearlab,squeeze(prctile(sig_post(:,ii-1,:),[50],3)),'color',[.7  0  0],'LineWidth',1.5);
        alpha(.5);
        legend(ggg,sprintf('Median'),sprintf('68%% Posterior'))
        title('Legend','FontSize',FontSize,'FontWeight','Bold')
        set(gca,'FontSize',FontSize,'FontWeight','Bold','Xlim',[0 1],'Ylim',[0 1],'XTickLabel',[],'YTickLabel',[],'Box','On')
    else
        tmpb = [squeeze(prctile(sig_post(:,ii,:),16,3));flipdim(squeeze(prctile(sig_post(:,ii,:),84,3)),1)];
        tmp_min=floor(min(tmpb(:))*10)/10;
        tmp_max=ceil(max(tmpb(:))*10)/10;
        patch(tmpa,tmpb,[255 182 193]./255,'EdgeColor',[255 182 193]./255)
        plot(yearlab,squeeze(prctile(sig_post(:,ii,:),[50],3)),'color',[.7  0  0],'LineWidth',1.0)
        axis('tight');grid;alpha(.5);
        set(gca,'Fontsize',Fontsize,'YLim',[tmp_min tmp_max],'Box','On');
        title(varname(ii),'FontSize',Fontsize,'FontWeight','Bold')
        xlabel(sprintf('\n'))
    end
end
saveTightFigure(gcf,sprintf('%sVolatility_Model_%d_Prior_%d_Level_2BR.tiff',save_path,data_sel,prior_sel))
delete(sprintf('%sVolatility_Model_%d_Prior_%d_Level_2BR.tiff',save_path,data_sel,prior_sel))
print(gcf,'-dtiff','-r1200',sprintf('%sVolatility_Model_%d_Prior_%d_Level_2BR',save_path,data_sel,prior_sel))
close('all')

%==========================================================================
% Forecast Moments calculation & figure creation
%==========================================================================
number_of_draws = nrep ;
for horz = 20 % forecast horizon: 5 year ahaed
    
    forecast.mean   = zeros(M,t,number_of_draws);
    forecast.mse    = zeros(M,M,t,number_of_draws);
    forecast.corr   = zeros(M,M,t,number_of_draws);
    y_post          = repmat(y_obs',[1 1 nrep])-me_post_diff;     % Data ordering: Y PIE R_Short M2 M0 R_Long
    
    tic
    % -------------------------------------------------------------------------
    % Posterior simulation of forecast moments
    % -------------------------------------------------------------------------
    for l = 1 : number_of_draws
        if mod(l,10)==0; fprintf('Draw:\t%6.0f\t(%6.2f min)\n',l,toc/60);end
        
        ylag    = mlag2(y_post(:,:,l),p);% y_post is [T x M x nrep]. ylag is [T x (Mp)]
        ctemp1  = zeros(M,M*p);
        Btdraw  = Bt_post(:,:,l);
        capAt   = zeros(M*t,M);
        for i = 1 : t
            capatemp    = eye(M);
            aatemp      = At_post(:,i,l);
            ic = 1;
            for j = 2 : M
                capatemp(j,1:j-1) = aatemp(ic:ic+j-2,1)';
                ic = ic + j - 1;
            end
            capAt((i-1)*M+1:i*M,:) = capatemp;
        end
        sigtemp = eye(M);
        sigt    = zeros(M*t,M);
        for i = 1:t
            for j = 1:M
                sigtemp(j,j) = exp(0.5*Sigt_post(j,i,l));
            end
            sigt((i-1)*M+1:i*M,:) = sigtemp;
        end
        Ht      = zeros(M*t,M);
        Htsd    = zeros(M*t,M);
        counter = 0;
        J       = [eye(M) zeros(M,M*(p-1))];
        for i = 1:t;
            stem        = sSiigt((i-1)*M+1:i*M,:);
            Hsd         = capAt((i-1)*M+1:i*M,:)\stem;
            Hdraw       = Hsd*Hsd';
            intercepts  = Btdraw(1:M,i);
            A1          = reshape(Btdraw(M+1:M+M^2,i),M,M)';
            A2          = reshape(Btdraw(M+M^2+1:end,i),M,M)';
            ctemp1      = [A1 A2; eye(M*(p-1)) zeros(M*(p-1),M)];
            tmp         = ylag(i,:)';
            
            tmp1=0; tmp2=0;
            for k=1:horz
                tmp     = [intercepts;zeros(M,1)]+ctemp1*tmp;
                phi     = J*ctemp1^(horz-1)*J';
                tmp1    = tmp1+phi*Hdraw*phi';
                tmp2    = tmp2+(phi*Hsd).^2;
            end;
            % horizon horz posterior forecast moments
            forecast.mean(:,i,l)=tmp(1:5,1);  % forecast mean vector at time i of draw l (VAR ordering: y, pie, r, sp, m)
            % dimension (M x T x nrep)
            forecast.mse(:,:,i,l)=tmp1; % forecast MSE matrix at time i of draw l (VAR ordering)
            % dimension (M x M x T x nrep)
            S=diag(diag(forecast.mse(:,:,i,l)));
            forecast.corr(:,:,i,l)=(S^(-1/2))*forecast.mse(:,:,i,l)*(S^(-1/2)); % forecast correlation matrix
            % dimension (M x M x T x nrep)
        end
    end
    
    forecast.mean   = prctile(forecast.mean,[5 10 16 20:5:80 84 90 95],3);
    forecast.mse    = prctile(sqrt(forecast.mse) ,[5 10 16 20:5:80 84 90 95],4);
    forecast.corr   = prctile(forecast.corr,[5 10 16 20:5:80 84 90 95],4);
    for  ii = 1 : M
        forecast.std(ii,:,:) = forecast.mse(ii,ii,:,:);
    end
    
    clear forecast.mse
    
    eval(sprintf('save %sCondMoments_AAMW_ME_Model_%d_Prior_%d_Horizon_%d forecast ',save_path,data_sel,prior_sel,horz));
    
    % -------------------------------------------------------------------------
    % Figure 3: (Forecast-Based)Conditional Mean
    % -------------------------------------------------------------------------
    FontSize = 4;
    tmp_axis = [yearlab;flipdim(yearlab,1)];
    RGB_Color = [173 216 230;
        0 191 255]./255;
    clear tmp_handle g
    figure('Name','Forecast Means');% orient('Landscape')
    for ii = 1 : M+1
        tmp_handle(ii) = subplot(3,2,ii);hold('on')
        if ii>M
            jjj=1;
            for jj = 3;%[1 3];%1 : 9
                g(1)=patch(tmp_axis,tmp_series,RGB_Color(jjj,:),'EdgeColor',RGB_Color(jjj,:));
                jjj=jjj+1;
            end
            g(2)=plot(yearlab,squeeze(forecast.mean(ii-1,:,10)),'Color',[65 105 225]/255,'LineWidth',2.0);
            alpha(.5);set(gca,'XLim',[0 1],'XTickLabel',[],'YLim',[0 1],'YTickLabel',[])
            title(sprintf('Legend'),'FontSize',FontSize,'FontWeight','Bold')
            legend(g,sprintf('68%% Posterior'),sprintf('Posterior Median'),'Location','NorthEast')
        else
            jjj=1;
            for jj = 3;%[1 3];%1 : 9
                tmp_series = [squeeze(forecast.mean(ii,:,jj))';flipdim(squeeze(forecast.mean(ii,:,end+1-jj))',1)];
                % patch(tmp_axis,tmp_series,RGB_Color(jjj,:),'EdgeColor',RGB_Color(jjj+1,:))
                patch(tmp_axis,tmp_series,RGB_Color(jjj,:),'EdgeColor',RGB_Color(jjj,:))
                jjj=jjj+1;
            end
            plot(yearlab,squeeze(forecast.mean(ii,:,10))*0,'Color',[0 0 0]/255,'LineWidth',1.0)
            plot(yearlab,squeeze(forecast.mean(ii,:,10)),'Color',[65 105 225]/255,'LineWidth',1.0)
            grid('on');axis('tight');alpha(.5)
            title(sprintf('%s',char(varname(ii))),'FontSize',FontSize,'FontWeight','Bold')
            xlabel(sprintf('\n'))
        end
    end
    set(tmp_handle,'FontSize',FontSize,'LineWidth',.5,'Box','On')
    
    saveTightFigure(gcf,sprintf('%sConditionalMean_Horizon_%d_Model_%d_Prior_%d.tiff',save_path,horz,data_sel,prior_sel))
    delete(sprintf('%sConditionalMean_Horizon_%d_Model_%d_Prior_%d.tiff',save_path,horz,data_sel,prior_sel))
    print(gcf,'-dtiff','-r1200',sprintf('%sConditionalMean_Horizon_%d_Model_%d_Prior_%d',save_path,horz,data_sel,prior_sel))
    

    % -------------------------------------------------------------------------
    % Figure 4: (Forecast-Based) Correlation
    % -------------------------------------------------------------------------
    tmp_1 = reshape(1:M^2,M,M)';
    tmp_2 = tril(reshape(1:M^2,M,M)',-1);
    tmp_sel = nonzeros(tril(reshape(1:M^2,M,M)',-1));
    
    clear tmp_handle g
    figure('Name','Forecast Correlation');% orient('Landscape')
    for ii = 1 : length(tmp_sel)+1
        if ii>length(tmp_sel)
            tmp_handle(ii) = subplot(3,4,[ii:ii+1]);hold('on');
            jjj=1;
            for jj = 3;%[1 3];%1 : 9
                g(1)=patch(tmp_axis,tmp_series,RGB_Color(jjj,:),'EdgeColor',RGB_Color(jjj+1,:));
            end
            g(2)=plot(yearlab,squeeze(forecast.corr(row_i,col_i,:,10)),'Color',[65 105 225]/255,'LineWidth',2.0);
            alpha(.5);set(gca,'XLim',[0 1],'XTickLabel',[],'YLim',[0 1],'YTickLabel',[])
            title(sprintf('Legend'),'FontSize',FontSize,'FontWeight','Bold')
            legend(g,sprintf('68%% Posterior'),sprintf('Posterior Median'),'Location','NorthEast')
        else
            tmp_handle(ii) = subplot(3,4,ii);hold('on');
            jjj=1;
            [row_i,col_i]=find(tmp_2==tmp_sel(ii));
            for jj = 3;%[1 3];%1 : 9
                tmp_series = [squeeze(forecast.corr(row_i,col_i,:,jj));flipdim(squeeze(forecast.corr(row_i,col_i,:,end+1-jj)),1)];
                patch(tmp_axis,tmp_series,RGB_Color(jjj,:),'EdgeColor',RGB_Color(jjj+1,:))
                jjj=jjj+1;
            end
            plot(yearlab,0*squeeze(forecast.corr(row_i,col_i,:,10)),'Color',[0 0 0],'LineWidth',1.0)
            plot(yearlab,squeeze(forecast.corr(row_i,col_i,:,10)),'Color',[65 105 225]/255,'LineWidth',1.0)
            grid('on');axis('tight');alpha(.5)
            title(sprintf('%s vs. %s',char(varname(row_i)),char(varname(col_i))),'FontSize',FontSize,'FontWeight','Bold')
            xlabel(sprintf('\n'))
        end
    end
    set(tmp_handle,'FontSize',FontSize,'LineWidth',.5,'YLim',[-1 1],'Box','On')
    
    saveTightFigure(gcf,sprintf('%sConditionalCorr_Horizon_%d_Model_%d_Prior_%d.tiff',save_path,horz,data_sel,prior_sel))
    delete(sprintf('%sConditionalCorr_Horizon_%d_Model_%d_Prior_%d.tiff',save_path,horz,data_sel,prior_sel))
    print(gcf,'-dtiff','-r1200',sprintf('%sConditionalCorr_Horizon_%d_Model_%d_Prior_%d',save_path,horz,data_sel,prior_sel))
    
    close all
    
end

close all;
