%% ========================================================================
% Main Code to run a tvp-var featuring measurement errors in the 
% obvervables developed for the paper: 
% 
%       "Drifts and Volatilities under Measurement Error:
%        Assessing Monetary Policy Shocks over the Last Century"
% 
%   to appear in Quantitative Economics
% 
%   by 
%   Pooyan Amir-Ahmadi (amir@econ.uni-frankfurt.de)
%   Christian Matthes (christian.matthes@gmail.com)
%   Mu-Chun Wang (Mu-Chun.Wang@wiso.uni-hamburg.de)
% 
% -------------------------------------------------------------------------
% Note: Codes are based on codes provided by Korrobilis and Koop for
%       tvp-var. The read me file contains all information on modication 
%       and extensions. 
% 
% If you wish to adjust this code for your data and application also search
% do case-sensitive sarch or "NOTE" to track which part is hard coded and 
% therefore needs modification and adjustment to user data. Note that this 
% won't be an exhausitve list. 
% 
%==========================================================================

%% 0. Housekeeping
% =========================================================================
clear('all'); close('all'); clc;

randn('state', sum(100*clock));
rand('twister',sum(100*clock));

main_path = [pwd '\'];
save_path = [pwd '\Results\'];

addpath([pwd '\Data\'])
addpath([pwd '\MiscCodes\'])

cd(main_path)


%% I. Load data and set specifications
% =========================================================================

load newnewdata_1876Q1_2011Q2_yoy   % Data ordering: Y PIE R_Short M2 M0 R_Long
                                    % [Real GDP Growth, Inflation, short rate, 
                                    % M2 growth, Money grwoth, long rate ]

spread      = data(:,6)-data(:,3);  
data_sel    = 2;                    % 1: M2, 2: M0; default is M0 growth, alternatively you can run the model based on M2 growth
prior_sel   = 1;                    % 1: k_q = .01, 2: k_q = .05, 3: k_q = .10, 4: k_q = sqrt(3.5e-4);
                                    % Here you can choose the tuning parameter fo the extend 

switch data_sel
    case 1 % VAR ordering: Y PIE R_Short Spread M2
        Y = [data(:,1) data(:,2) data(:,3) spread data(:,4)];
    case 2 % VAR ordering: Y PIE R_Short Spread M0
        Y = [data(:,1) data(:,2) data(:,3) spread data(:,5)];
end

t           = size(Y,1);            % t is the time-series observations of Y
M           = size(Y,2);            % M is the dimensionality of Y
tau         = 152;                  % size of the training sample 1876:1 - 1913:4
p           = 2;                    % # lags
numa        = M*(M-1)/2;            % non-0,non-1 elements of A_t

% VAR EQUATION
% -------------------------------------------------------------------------
ylag        = mlag2(Y,p);           % lagged Y is [T x M]. ylag is [T x (Mp)]
ylag        = ylag(p+tau+1:t,:);    % Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T
K           = M + p*(M^2);          % K is the number of elements in the state vector
Z           = zeros((t-tau-p)*M,K); % Create Z_t matrix.
for i = 1:t-tau-p
    ztemp = eye(M);
    for j = 1:p
        xtemp = ylag(i,(j-1)*M+1:j*M);
        xtemp = kron(eye(M),xtemp);
        ztemp = [ztemp xtemp];
    end
    Z((i-1)*M+1:i*M,:) = ztemp;
end
y           = Y(tau+p+1:t,:)';
y_obs       = y;                    % Note: this is OBSERVED y
yearlab     = yearlab(tau+p+1:end);
t           = size(y,2);            % t is now smaller/adjusted

% Specify posterior sampler
% -------------------------------------------------------------------------
nrep        =  2000;    % Number of draws kept from stationary phase
nburn       = 98000;    % Number of burn-in draws
nsave       =  2000;    % Number of draws saved in each file
nskip       =     1;    % Number of draws skipped
it_print    =   100;    % Print in command window every "it_print"-th iteration

draw_index  = randsample(nrep,nsave)+nburn; % If you wish to select a random subset of draws for further evaluation 
% draw_index  = 1:nsave;                    

DO_TRAIN    = 1;    % for initialization set 1: to use training sample, 2: to use uninformative prior
DO_SAVE     = 1;    % Save figures and results
DO_OLDVAL   = 0;	% 1: start at last draw of previous run; 0: otherwise
DO_STAB     = 0;    % 1: check for local stationarity, 0: otherwise (default)

% For estimation of measurement error process in the observables 
DO_PRIOR    = 2;    % 1: N-iG prior; 2: independent prior (our default choice)

% PRIORS
% -------------------------------------------------------------------------
switch DO_TRAIN
    case 1 % Training sample prior
        [B_OLS,VB_OLS,A_OLS,sigma_OLS,VA_OLS]= ts_prior(Y,tau,M,p);
    case 2 % Uninformative values
        A_OLS       = zeros(numa,1);
        B_OLS       = zeros(K,1);
        VA_OLS      = eye(numa);
        VB_OLS      = eye(K);
        sigma_OLS   = 0*ones(M,1);
end

switch prior_sel % Degree of varation in coefficients
    case 1 % Primiceri
        k_Q = 0.01;
    case 2 % Cogley & Sargent specification
        k_Q = sqrt(3.5e-4);
    case 3 % user definition
        % k_Q = ...
end
k_S     = 0.1;              % Degree of varation in covariance states
k_W     = 0.01;             % Degree of varation in SV

sizeW   = M;                % Size of matrix W
sizeS   = 1:M;              % Size of matrix S
sizeB   = K;                % Size of matrxi B


% Parameterization of Measurement Error Process
% -------------------------------------------------------------------------
y_train_stds    = std(Y(1:tau,:))';     % mean observed variables in training sample
y_train_means   = mean(Y(1:tau,:))';    % std of observed variables in training sample

% Prior for measurement error parameters: rhome and sigme2
% -------------------------------------------------------------------------
rhome_prmean    = zeros(M,1);
rhome_prvar     = ones(M,1)*(0.15^2);

sigme2_mode     = (0.25*y_train_stds).^2;       % mode of sigme2 (25% of the std of the training sample)
sigme2_a        = ones(M,1)*2;                  % inverse gamma degree of freedom, this is equal to the shape parameter a
sigme2_b        = sigme2_mode.*(sigme2_a+1);    % Here we reparameterze the inverse gamma mode into standard inverse gamma
                                                % specification Inv-Gamma(a,b)

% Kalman filter initial conditions for B(t), A(t) and log(SIGMA(t)).
B_0_prmean      = B_OLS;                % B_0 ~ N(B_OLS, 4Var(B_OLS))
B_0_prvar       = 4*VB_OLS;

A_0_prmean      = A_OLS;                % A_0 ~ N(A_OLS, 4Var(A_OLS))
A_0_prvar       = 4*VA_OLS;

sigma_prmean    = sigma_OLS;            % log(sigma_0) ~ N(log(sigma_OLS),I_n)
sigma_prvar     = 4*eye(M);

Q_prmean        = ((k_Q)^2)*tau*VB_OLS;
Q_prvar         = tau;

W_prmean        = ((k_W)^2)*(1 + sizeW)*eye(M); % W_prmean = 0.00001*eye(M);
W_prvar         = 1 + sizeW;                    % W_prvar = tau;

S_prmean        = cell(M-1,1);
S_prvar         = zeros(M-1,1);

ind = 1;
for ii = 2:M % S is block diagonal as in Primiceri (2005)
    S_prmean{ii-1}  = ((k_S)^2)*(1 + sizeS(ii-1))*VA_OLS(((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind);
    S_prvar(ii-1)   = 1 + sizeS(ii-1);
    ind             = ind + ii;
end

% Parameters of the 7 component mixture approximation to a log(chi^2) density:
q_s     = [  0.00730;  0.10556;  0.00002; 0.04395; 0.34001; 0.24566;  0.25750]; % probabilities
m_s     = [-10.12999; -3.97281; -8.56686; 2.77786; 0.61942; 1.79518; -1.08819]; % means
u2_s    = [  5.79596;  2.61369;  5.17950; 0.16735; 0.64009; 0.34023;  1.26261]; % variances

% INITIALIZE MATRICES:
% -------------------------------------------------------------------------
% Specify covariance matrices for measurement and state equations
consQ   = 0.0001;
consS   = 0.0001;
consH   = 0.01;
consW   = 0.0001;

Ht          = kron(ones(t,1),consH*eye(M));         % Initialize Htdraw, a draw from the VAR covariance matrix
Htchol      = kron(ones(t,1),sqrt(consH)*eye(M));   % Cholesky of Htdraw defined above
Qdraw       = consQ*eye(K);                         % Initialize Qdraw, a draw from the covariance matrix Q
Sdraw       = consS*eye(numa);                      % Initialize Sdraw, a draw from the covariance matrix S
Sblockdraw  = cell(M-1,1);                          % ...and then get the blocks of this matrix (see Primiceri)
ijc = 1;
for jj=2:M
    Sblockdraw{jj-1} = Sdraw(((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc);
    ijc = ijc + jj;
end
Wdraw       = consW*eye(M);                 % Initialize Wdraw, a draw from the covariance matrix W
Btdraw      = zeros(K,t);                   % Initialize Btdraw, a draw of the mean VAR coefficients, B(t)
Atdraw      = zeros(numa,t);                % Initialize Atdraw, a draw of the non 0 or 1 elements of A(t)
Sigtdraw    = zeros(M,t);                   % Initialize Sigtdraw, a draw of the log-diagonal of SIGMA(t)
sigt        = kron(ones(t,1),0.01*eye(M));  % Matrix of the exponent of Sigtdraws (SIGMA(t))
statedraw   = 5*ones(t,M);                  % initialize the draw of the indicator variable (of 7-component mixture of Normals approximation)
Zs          = kron(ones(t,1),eye(M));
prw         = zeros(numel(q_s),1);

rhomedraw   = zeros(M,t);                   % Initialize persistence parameters of measurement error shocks
sigme2draw  = sqrt(eps)*ones(M,t);          % Initialize variances of measurement error shocks
medraw      = zeros(t,M);                   % Initialize measurement error shocks, drawing from priors modes

% me_state_MSE = diag(sigme2_mode);         % Initialize initial guess measurement error MSE (case of measurement error 
                                            % directly on observables in growth rates)

me_state_MSE = blkdiag(eye(5)*sigme2_mode(1,1),eye(5)*sigme2_mode(2,1),sigme2_mode(3,1),...
             sigme2_mode(4,1),eye(5)*sigme2_mode(5,1)); % Initialize initial measurement error MSE

% Memory allocation
% -------------------------------------------------------------------------
Bt_post     = zeros(K,t,nsave);         % regression coefficients B(t)
At_post     = zeros(numa,t,nsave);      % lower triangular matrix A(t)
Sigt_post   = zeros(M,t,nsave);         % diagonal std matrix SIGMA(t)
Q_post      = zeros(K,K,nsave);         % covariance matrix Q of B(t)
ikc = 1;
for kk = 2:M
    Sdraw(((kk-1)+(kk-3)*(kk-2)/2):ikc,((kk-1)+(kk-3)*(kk-2)/2):ikc)=Sblockdraw{kk-1};
    ikc = ikc + kk;
end
S_post      = zeros(numa,numa,nsave);   % covariance matrix Q of B(t)
W_post      = zeros(M,M,nsave);         % covariance matrix Q of B(t)
sig_post    = zeros(t,M,nsave);         % diagonal of the VAR covariance matrix
cor_post    = zeros(t,numa,nsave);      % off-diagonal elements of the VAR cov matrix
rhome_post  = zeros(M,t,nsave);         % posteiror of persistence of measurement error process 
sigme2_post = zeros(M,t,nsave);         % posteiror of variance of measurement error process 
me_post     = zeros(t,M,nsave);         % posteiror of measurement error

% -------------------------------------------------------------------------
% Start from last draw of previous run if selected
% -------------------------------------------------------------------------
if DO_OLDVAL
    try
        eval(sprintf('load %sOldInitialization_AAMW_ME_Model_%d_Prior_%d',save_path,data_sel,prior_sel));
        disp('SUCESS')
        pause(1)
    catch
        disp('===========================================================')
        warning('Could not find starting values from previous run!')
        disp('===========================================================')
        pause(1)
    end
end


%% II. Start posterior sampler
%==========================================================================
tic; 
disp('Starting estimation...');
iirep       = 1;

for irep = 1 : nrep + nburn    % GIBBS iterations starts here
    if mod(irep,it_print) == 0;
        fprintf('Draw:\t%8.0f (took %6.2f minutes)\n',irep,toc/60);
    end
    
    % ---------------------------------------------------------------------
    % p(B_t,Q|Y,...): sample coeff. states and respective residual covariance
    % ---------------------------------------------------------------------
    draw_beta_corrected
    
    % ---------------------------------------------------------------------
    % p(A_t,S|Y,...): sample cov. states and respective residual covariance
    % ---------------------------------------------------------------------
    draw_alpha_corrected
    
    % ---------------------------------------------------------------------
    % Note: In the following step the crrigendum following DelNegro and Primiceri (2015, REStud) implemented
    %       More recent version of codes by Koop and Korrobilis have corrected for this!
    % ---------------------------------------------------------------------
    [statedraw,yss,capAt]           = draw_sTcomp_corrected(Atdraw,Sigtdraw,yhat,m_s,u2_s,q_s,M,t);
    [Sigtdraw,log_lik3,sigt,Wdraw]  = draw_sigma_correctedSep(statedraw,Wdraw,yss,Zs,m_s,u2_s,M,t,sigma_prmean,sigma_prvar,W_prmean,W_prvar);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dim_state   = M*p+5+5+1+1+5;        % dimension of the state vector = size of companion matrix  + # measurement errors (with lags)
                                        % NOTE: This is hard coded and needs to be adjusted different data set 
    Ht          = zeros(M*t,M);
    Htsd        = zeros(M*t,M);
    Qt          = zeros(dim_state,dim_state,t);
    for i = 1 : t
        inva                    = inv(capAt((i-1)*M+1:i*M,:));
        stem                    = sigt((i-1)*M+1:i*M,:);
        Hsd                     = capAt((i-1)*M+1:i*M,:)\stem;
        Hdraw                   = Hsd*Hsd';
        Ht((i-1)*M+1:i*M,:)     = Hdraw;  % H(t),
        Htsd((i-1)*M+1:i*M,:)   = Hsd;  % Cholesky of H(t)
        
        % forming state innovation covariance matrix Qt 
        % NOTE: This is hard coded and needs to be adjusted different data set 
        Qt(:,:,i)=blkdiag(Hdraw,zeros(M*(p-1)),diag([sigme2draw(1,i) zeros(1,4)]),...
            diag([sigme2draw(2,i) zeros(1,4)]),sigme2draw(3,i),sigme2draw(4,i),diag([sigme2draw(5,i) zeros(1,4)]));
    end
    
    % forming VAR(1) companion matrix
    biga = zeros(M*p,M*p);
    for j = 1:p-1
        biga(j*M+1:M*(j+1),M*(j-1)+1:j*M) = eye(M);
    end
    
    Ft          = zeros(dim_state,dim_state,t);
    intercepts  = zeros(dim_state,t);
    for i = 1:t %Get impulses recurssively for each time period
            % NOTE: This is hard coded and needs to be adjusted different data set
            intercepts(:,i) = [Btdraw(1:M,i);zeros(M*(p-1),1);zeros(17,1)]; % 
            bbtemp          = Btdraw(M+1:K,i);  % get the draw of B(t) at time i=1,...,T  (exclude intercept)
            splace          = 0;
            for ii = 1:p
                for iii = 1:M
                    biga(iii,(ii-1)*M+1:ii*M) = bbtemp(splace+1:splace+M,1)';
                    splace = splace + M;
                end
            end
             % forming state law of motion matrix Ft (Mxp+M,Mxp+M,t)
             % NOTE: This is hard coded and needs to be adjusted different data set
             rho1 = [rhomedraw(1,i) zeros(1,3); eye(4)];
             rho1 = [rho1 zeros(5,1)];
             rho2 = [rhomedraw(2,i) zeros(1,3); eye(4)];
             rho2 = [rho2 zeros(5,1)];
             rho5 = [rhomedraw(5,i) zeros(1,3); eye(4)];
             rho5 = [rho5 zeros(5,1)];
             
             Ft(:,:,i) = blkdiag(biga, rho1, rho2,rhomedraw(3,i),rhomedraw(4,i),rho5);
    end
    % NOTE: This is hard coded and needs to be adjusted different data set
    % forming selection matrix for the measurement equation
    J=[eye(M) zeros(M*(p-1)) zeros(M,17)]; 
    % select me for gdp growth
    J(1,M*p+1)=1;J(1,M*p+5)=-1;
    % select me for inflation
    J(2,M*p+6)=1;J(2,M*p+10)=-1;
    % select me for interest rate
    J(3,M*p+11)=1;
    % select me for spread
    J(4,M*p+12)=1;
    % select me for money growth
    J(5,M*p+13)=1;J(5,M*p+17)=-1;
    
   
    % NOTE: This is hard coded and needs to be adjusted different data set
    % initializing states with training sample variances of observables and unconditional variances of measurement errors
    state_0     = [Y(tau+1,:)';Y(tau,:)';zeros(17,1)];
    state_MSE_0 = blkdiag(diag(repmat(y_train_stds.^2,2,1)),me_state_MSE); 
    
    R           = zeros(M);
    
    % drawing unobserved TRUE y
    state_draw  = carter_kohn_me_with_const(y_obs,Qt,Ft,J,R,intercepts,dim_state,t,state_0,state_MSE_0);                        
    
    % picking up the correct ME from state vector
    % NOTE: This is hard coded and needs to be adjusted different data set
    medraw      = [state_draw(M*p+1,:)' state_draw(M*p+6,:)' state_draw(M*p+11,:)' state_draw(M*p+12,:)' state_draw(M*p+13,:)']; 
    
    % update smoothed unobserved data vector
    y           = state_draw(1:M,:); 
    
    % building Z regressor matrix, filling first p-lag inital values with
    % observed data, (assuming initial measurement errors = 0)
    Ys      = [Y(1:tau+p,:); y'];
    ylag    = mlag2(Ys,p);                      % lagged Y is [T x M]. ylag is [T x (Mp)]
    ylag    = ylag(p+tau+1:end,:);              % Form RHS matrix X_t = [1 y_t-1 y_t-2 ... y_t-k] for t=1:T
    Z       = zeros((size(Ys,1)-tau-p)*M,K);    % Create Z_t matrix.
    for i = 1 : size(Ys,1)-tau-p
        ztemp = eye(M);
        for j = 1:p
            xtemp = ylag(i,(j-1)*M+1:j*M);
            xtemp = kron(eye(M),xtemp);
            ztemp = [ztemp xtemp];
        end
        Z((i-1)*M+1:i*M,:) = ztemp;
    end
    
    lag_medraw =lag(medraw,1);
    

    %% Measurement error sampling
    % =====================================================================
    % NOTE: This is hard coded and needs to be adjusted different data set
    
    % ---------------------------------------------------------------------
    % first variable Y
    % ---------------------------------------------------------------------
    % 1st Break date
    break_date1 = find(yearlab==1930); % until 1917Q4
    
    [rhomedraw_tmp, sigme2draw_tmp]=...
        norm_indp_invgamm(medraw(2:break_date1,5),lag_medraw(2:break_date1,1),rhome_prmean(1,1),rhome_prvar(1,1),sigme2_a(1,1),sigme2_b(1,1),rhomedraw(1,1));
    
    sigme2draw(1,1:break_date1) = sigme2draw_tmp;
    rhomedraw(1,1:break_date1)  = rhomedraw_tmp;
    
    % 2nd Break date
    break_date2 = 130; % until 1935Q4
    
    [rhomedraw_tmp, sigme2draw_tmp]=...
        norm_indp_invgamm(medraw(break_date1+1:break_date2,1),lag_medraw(break_date1+1:break_date2,1),rhome_prmean(1,1),rhome_prvar(1,1),sigme2_a(1,1),sigme2_b(1,1),rhomedraw(1,break_date1+1));
    
    sigme2draw(1,break_date1+1:break_date2) = sigme2draw_tmp;
    rhomedraw(1,break_date1+1:break_date2)  = rhomedraw_tmp;
    
    
    % ---------------------------------------------------------------------
    % 2nd variable PIE
    % ---------------------------------------------------------------------
    % 1st Break date
    break_date = 130; % until 1946Q4
    [rhomedraw_tmp, sigme2draw_tmp]=...
        norm_indp_invgamm(medraw(2:break_date,2),lag_medraw(2:break_date,2),rhome_prmean(2,1),rhome_prvar(2,1),sigme2_a(2,1),sigme2_b(2,1),rhomedraw(2,1));
    
    sigme2draw(2,1:break_date)  = sigme2draw_tmp;
    rhomedraw( 2,1:break_date)  = rhomedraw_tmp;
    sigme2draw(2,break_date:end)= sqrt(eps);
    rhomedraw( 2,break_date:end)= 0;
    
    
    % ---------------------------------------------------------------------
    % 3rd variable R
    % ---------------------------------------------------------------------
    % 1st Break date
    break_date = 22; % until 1919Q4
    
    [rhomedraw_tmp, sigme2draw_tmp]=...
        norm_indp_invgamm(medraw(2:break_date,3),lag_medraw(2:break_date,3),rhome_prmean(3,1),rhome_prvar(3,1),sigme2_a(3,1),sigme2_b(3,1),rhomedraw(3,1));
    
    sigme2draw(3,1:break_date)  = sigme2draw_tmp;
    rhomedraw( 3,1:break_date)  = rhomedraw_tmp;
    sigme2draw(3,break_date:end)= sqrt(eps);
    rhomedraw( 3,break_date:end)= 0;
    
    
    % ---------------------------------------------------------------------
    % 4th variable Spread
    % ---------------------------------------------------------------------
    % 1st Break date
    break_date = 22; % until 1919Q4
    
    [rhomedraw_tmp, sigme2draw_tmp]=...
        norm_indp_invgamm(medraw(2:break_date,4),lag_medraw(2:break_date,4),rhome_prmean(4,1),rhome_prvar(4,1),sigme2_a(4,1),sigme2_b(4,1),rhomedraw(4,1));
    
    sigme2draw(4,1:break_date)  = sigme2draw_tmp;
    rhomedraw( 4,1:break_date)  = rhomedraw_tmp;
    sigme2draw(4,break_date:end)= sqrt(eps);
    rhomedraw( 4,break_date:end)= 0;
    
    
    % ---------------------------------------------------------------------
    % 5th variable M0
    % ---------------------------------------------------------------------
    % 1st Break date
    break_date1 = 14; % until 1917Q4
    
    [rhomedraw_tmp, sigme2draw_tmp]=...
        norm_indp_invgamm(medraw(2:break_date1,5),lag_medraw(2:break_date1,5),rhome_prmean(5,1),rhome_prvar(5,1),sigme2_a(5,1),sigme2_b(5,1),rhomedraw(5,1));
    
    sigme2draw(5,1:break_date1) = sigme2draw_tmp;
    rhomedraw( 5,1:break_date1) = rhomedraw_tmp;
    
    % 2nd Break date
    break_date2 = 86; % until 1935Q4
    
    [rhomedraw_tmp, sigme2draw_tmp]=...
        norm_indp_invgamm(medraw(break_date1+1:break_date2,5),lag_medraw(break_date1+1:break_date2,5),rhome_prmean(5,1),rhome_prvar(5,1),sigme2_a(5,1),sigme2_b(5,1),rhomedraw(5,break_date1+1));
    
    sigme2draw(5,break_date1+1:break_date2) = sigme2draw_tmp;
    rhomedraw( 5,break_date1+1:break_date2) = rhomedraw_tmp;
    
    % 2nd Break date
    break_date3 = 178; % until 1958Q4
    
    [rhomedraw_tmp, sigme2draw_tmp]=...
        norm_indp_invgamm(medraw(break_date2+1:break_date3,5),lag_medraw(break_date2+1:break_date3,5),rhome_prmean(5,1),rhome_prvar(5,1),sigme2_a(5,1),sigme2_b(5,1),rhomedraw(5,break_date2+1));
    
    sigme2draw(5,break_date2+1:break_date3) = sigme2draw_tmp;
    rhomedraw( 5,break_date2+1:break_date3) = rhomedraw_tmp;
    
    sigme2draw(5,break_date3:end)=sqrt(eps);
    rhomedraw(5,break_date3:end)=0;

    
    % ---------------------------------------------------------------------
    % Start saving after burn-in
    % ---------------------------------------------------------------------
    % if irep > nburn && mod(irep,nskip)==0; 
    if ismember(irep,draw_index) == 1;
        % Get time-varying correlations and variances
        stemp6 = zeros(M,1);
        stemp5 = [];
        stemp7 = [];
        for i = 1:t
            stemp8 = corrvc(Ht((i-1)*M+1:i*M,:));
            stemp7a = [];
            ic = 1;
            for j = 1:M
                if j>1;
                    stemp7a = [stemp7a ; stemp8(j,1:ic)']; %#ok<AGROW>
                    ic = ic+1;
                end
                stemp6(j,1) = sqrt(Ht((i-1)*M+j,j));
            end
            stemp5 = [stemp5 ; stemp6'];
            stemp7 = [stemp7 ; stemp7a'];
        end
        sig_post(:,:,iirep)     = stemp5;   % diagonal of the VAR covariance matrix
        cor_post(:,:,iirep)     = stemp7;   % off-diagonal elements of the VAR cov matrix
        Bt_post(:,:,iirep)      = Btdraw;   % regression coefficients B(t)
        At_post(:,:,iirep)      = Atdraw;   % lower triangular matrix A(t)
        Sigt_post(:,:,iirep)    = Sigtdraw; % diagonal std matrix SIGMA(t)
        Q_post(:,:,iirep)       = Qdraw;    % covariance matrix Q of B(t)
        ikc = 1;
        for kk = 2:M
            Sdraw(((kk-1)+(kk-3)*(kk-2)/2):ikc,((kk-1)+(kk-3)*(kk-2)/2):ikc)=Sblockdraw{kk-1};
            ikc = ikc + kk;
        end
        S_post(:,:,iirep)       = Sdraw;    % covariance matrix S of A(t)
        W_post(:,:,iirep)       = Wdraw;    % covariance matrix W of SIGMA(t)
        me_post(:,:,iirep)      = medraw;   % sampled measurement error
        rhome_post(:,:,iirep)   = rhomedraw;
        sigme2_post(:,:,iirep)  = sigme2draw;
        
        iirep = iirep + 1;
    end % END saving after burn-in results
end %END main Gibbs loop (for irep = 1:nrep+nburn)

if DO_SAVE
	eval(sprintf('save %sOldInitialization_AAMW_ME_Level_Model_%d_Prior_%d Wdraw Sdraw Qdraw Sigtdraw Atdraw Btdraw sigt statedraw Sblockdraw medraw rhomedraw sigme2draw y Z Ht',save_path,data_sel,prior_sel));
    eval(sprintf('save %sAAMW_ME_Level_Model_%d_Prior_%d  ',save_path,data_sel,prior_sel));
end

fprintf('\n')
fprintf('====================================================================================\n')
fprintf('Running full code took:\t%8.2f\t hours!\n',toc/3600)
fprintf('====================================================================================\n')
fprintf('\n')



