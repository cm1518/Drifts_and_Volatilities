%% ========================================================================
% Main Code to run generilized impulse response functions in a  tvp-var 
% featuring measurement errors in the obvervables developed for the paper: 
% 
%       "Drifts and Volatilities under Measurement Error:
%        Assessing Monetary Policy Shocks over the Last Century"
% 
%   by 
%   Pooyan Amir-Ahmadi (amir@econ.uni-frankfurt.de)
%   Christian Matthes (christian.matthes@gmail.com)
%   Mu-Chun Wang (Mu-Chun.Wang@wiso.uni-hamburg.de)
% 
% -------------------------------------------------------------------------
% Note: Code requires to run Main_me_tvpvar.m first with mat files saved
%       in the folder <Results>.
% 
% #########################################################################
% #########################################################################
% 
% WARNING: RUNNING THIS SCRIPT CAN BE BOTH VERY TIME AND MEMORY INTENSIVE.
% 
% #########################################################################
% #########################################################################
% 
%==========================================================================

%% 0. Housekeeping
% =========================================================================
clear('all'); close('all'); clc;

% Make sure to have the same cd as the main script <Main_me_tvpvar.m>
main_path = [pwd '\'];
save_path = [pwd '\Results\'];

addpath(main_path)
addpath(save_path)
addpath([pwd '\Data\'])
addpath([pwd '\MiscCodes\'])

cd(main_path)

prior_sel   = 1;
data_sel    = 2;
result_path = save_path;

%% I. Load necessary Mat files
% =========================================================================
eval(sprintf('load %sAAMW_ME_Level_Model_%d_Prior_%d',result_path,data_sel,prior_sel));

number_of_draws = 250;
nalpha          = 50;                       % # of rotations per posterior draw
NMCInteg        = nalpha;
SC              = 0;                        % 0: No local stability conditions imposed (default)
final_kept_draw = number_of_draws*nalpha;

nhorz           = 21;                       % # GIRF horizon
h               = nhorz;
sr_horz         = 3;                        % # of horizon for sign restricitons to hold
DO_Y            = 1;                        % 0: agnostic w.r.t. output, 1: impose sr on output
sr_end          = zeros(number_of_draws,t);
draw_index      = randsample(nrep,number_of_draws);

tube.srgirf     = nan(M,t,nhorz,final_kept_draw);
tube.norm       = nan(M,t,nhorz,number_of_draws*3);

save_index      = [1:nalpha:final_kept_draw final_kept_draw+1];
save_index1     = [1:3:number_of_draws*3 number_of_draws*3+1];

tic;

for l = 1 : number_of_draws 
    fprintf('Draw: %8.0f\t(%4.2f\t minutes.)\n',l,toc/60)
    for i = 1 : t        
        if mod(i,1)==0;fprintf('Draw: %8.0f at time\t %3.0f\t(%4.2f\t minutes.)\n',l,i,toc/60);end

        g_rand = l;
        %% Here we do GIRFs á la Koop, Pesaran and Potter (1996)
        % =====================================================================
        Btmp    = Bt_post(:,i,g_rand);
        Atmp    = squeeze(At_post(:,i,g_rand));
        Sigtmp  = squeeze(Sigt_post(:,i,g_rand));
        Qtmp    = squeeze(Q_post(:,:,g_rand));
        Wtmp    = squeeze(W_post(:,:,g_rand));
        Stmp    = squeeze(S_post(:,:,g_rand));
        stem    = diag(sig_post(i,:,g_rand));
        nhor    = 21;
        tt      = p+1;          % for tt = p+1:length(yearlab)
        YYinit  = Y(length(Y)-t+1-p:length(Y),:)';
        if i <= p
            YY      = YYinit(:,i:i+p-1,:);
        else
            YY      = y(:,i-p:i-1);
        end
        [ind,girf,sirf]     = do_sr_girf_public(NMCInteg,SC,K,M,numa,nhor,p,Btmp,Atmp,Sigtmp,Qtmp,Wtmp,Stmp,YY,DO_Y,sr_horz);
        
        tmp1 = prctile( 0.25*(sirf ./ repmat(sirf(:,3,1),[1 M nhor])),[16 50 84],1); % norm impact interest rate response to .25 percent
        
        tube.srgirf(:,i,:,save_index(1):save_index(2)-1)    = permute(sirf,[2 3 1]);
        tube.norm(:,i,:,save_index1(1):save_index1(2)-1)    = permute(tmp1,[2 3 1]);
        
        
    end
    save_index(1)=[];
    save_index1(1)=[];
end

clear -regexp _post$;
clearvars -except tube data_sel prior_sel M t yearlab nhorz h sr_horz DO_Y sr_end draw_index nrep l

% Save full posterior SR identiffied GIRF 
% -------------------------------------------------------------------------
eval(sprintf('save AAMW_ME_Level_Model_%d_Prior_%d_Yres_%d_SRHorz_%d_sr_GIRF.mat',data_sel,prior_sel,DO_Y,sr_horz))

tube.srgirf     = squeeze(prctile(tube.srgirf,[5:5:95],4));
tube.norm       = squeeze(prctile(tube.norm,[5:5:95],4));

% Save specific percentiles of posterior SR identiffied GIRF 
% -------------------------------------------------------------------------
eval(sprintf('save 2AAMW_ME_Level_%d_Prior_%d_Yres_%d_SRHorz_%d_sr_GIRF.mat',data_sel,prior_sel,DO_Y,sr_horz))

fprintf('\n')
fprintf('==================================================================\n')
fprintf('Full GIRF Run: %8.0f\t(%4.2f\t minutes.)\n',l,toc/60)
fprintf('==================================================================\n')

