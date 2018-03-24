function [ind,girf,sirf] = do_sr_girf_public(NMCInteg,SC,K,M,numa,nhor,p,Btmp,Atmp,Sigtmp,Qtmp,Wtmp,Stmp,YY,DO_Y,sr_horz)


%% Part 1: Simulate TV-Parameters into the future:
% =========================================================================

% -------------------------------------------------------------------------
% Step 1: Simulate A_t, \Sigma_t and Ht
% -------------------------------------------------------------------------
for jj = 1 : nhor+p+1
    Atmp(:,:,jj+1)      = Atmp(:,:,jj) + (real(sqrtm(Stmp)) * randn(numa,1));       % Step 1: Simulate A_t
    Sigtmp(:,:,jj+1)    = Sigtmp(:,:,jj) + (Wtmp*randn(M,1));% om_e(:,jj);          % Step 2: Simulate \Sigma_t
    capatemp = eye(M);                                                              % Step 3: Simulate Q_t
    ic=1;
    for j = 2:M
        capatemp(j,1:j-1) = Atmp(ic:ic+j-2,:,jj)';
        ic = ic + j - 1;
    end
    inva            = inv(capatemp);
    Hsd             = inva*diag(Sigtmp(:,:,jj));
    VCVtmp(:,:,jj)  = Hsd*Hsd';
end

% -------------------------------------------------------------------------
% Step 1: Simulate B_t
% -------------------------------------------------------------------------
biga        = zeros(M*p,M*p);
for j = 1 : p-1
    biga(j*M+1:M*(j+1),M*(j-1)+1:j*M)   = eye(M);
end
B           = zeros(M,M*p+1,nhor+p+1);
jj          = 2; % Simulating the matrix B:
trial       = 0;
maxtrial    = 50;

%% Final Stuff
sirf = nan(NMCInteg,M,nhor);

while jj <= nhor+p+1
    bbtemp      = Btmp(:,:,jj-1) + real(mysqrt_varm(Qtmp))*randn(K,1);
    consttmp1   = bbtemp(1:M);%Btdrawc(M+1:K,i);  % get the draw of B(t) at time i=1,...,T  (exclude intercept)
    bbtemp(1:M) = [];
    splace  = 0;
    for ii = 1:p
        for iii = 1:M
            biga(iii,(ii-1)*M+1:ii*M)   = bbtemp(splace+1:splace+M,1)';
            splace = splace + M;
        end
    end
    consttmp        = zeros(numa,1);
    consttmp(1:M)   = consttmp1;
    b               = [consttmp biga];
    if SC==0
        Btmp(:,:,jj)=[consttmp1;bbtemp]; % sd(:,:,jj)=btemp+const;
        B(:,:,jj)=b(1:5,:);
        jj=jj+1;
    else
        if max(abs(eig(biga)))<1
            % % sd(:,:,jj)=dd; sd(:,:,jj)=btemp+const;
            Btmp(:,:,jj)=[consttmp1;bbtemp]; % sd(:,:,jj)=btemp+const;
            B(:,:,jj)=b;
            jj=jj+1;
            trial=0;
        elseif max(abs(eig(biga)))>1 && trial<maxtrial
            trial=trial+1;
        elseif max(abs(eig(biga)))>1 && trial==maxtrial
            sirf=[];
            dirf=[];
            aggirf=[];
            ind=0;
            d=0;
            return
        end
    end
end

%% Part 2: Calculate GIRFs
% =========================================================================
tel     = 1;
count   = 0;
maxcount= 2500;

while tel<NMCInteg && count<maxcount
    count = count+1;
    a0              = impact(VCVtmp(:,:,1),M);              % contemporaneous impact matrix
    Y_BEN1          = [YY zeros(M,nhor)];
    Y_MP            = Y_BEN1;   %
    SHOCKS          = randn(M,nhor);                        % same shocks for all identified shocks
    Y_BEN1(:,p+1)   = a0*SHOCKS(:,1);                       % The impact Matrix
    Y_MP(:,p+1)     = a0*[SHOCKS(1:2,1)' 1+SHOCKS(3,1) SHOCKS(4:end,1)']';
    
    for tt = p+2 : nhor+p
        shocks          = mysqrt_varm(VCVtmp(:,:,tt-p))*SHOCKS(:,tt-p);
        Y_BEN1(:,tt)    = B(:,:,tt-p)*[1; myvec(myfliplr(Y_BEN1(:,tt-p:tt-1)))]+shocks;
        Y_MP(:,tt)      = B(:,:,tt-p)*[1; myvec(myfliplr(Y_MP(:,tt-p:tt-1)))]  +shocks;
    end
    Y_BEN1      = Y_BEN1(:,p+1:size(Y_BEN1,2));   %3x20 matrices
    Y_MP        = Y_MP(:,p+1:size(Y_MP,2));
    irf_cand    = Y_MP-Y_BEN1;   % IRF after 1st shock: 3x20 matrix i.e. response of all 3 variables to first shock
    cirf_cand   = cumsum(irf_cand,2);
    
    
    %% 2: Check Sign restrictions here
    % ---------------------------------------------------------------------
    switch DO_Y
        case 0 % 1:'GDP', 2:'Inflation', 3:'Interest rate', 4:'Spread', 5:'Money growth'
            %% 4. Check Sign Restriction
            if all(cirf_cand(2,1:sr_horz)<=0) && all(irf_cand(3,1:sr_horz)>0) && all(irf_cand(5,1:sr_horz)<=0)
                SR_check = 1;
            elseif all(cirf_cand(2,1:sr_horz)>0) && all(irf_cand(3,1:sr_horz)<0) && all(irf_cand(5,1:sr_horz)>0)
                SR_check = -1;
            else
                SR_check = 0;
            end
        case 1% 1:'GDP', 2:'Inflation', 3:'Interest rate', 4:'Spread', 5:'M2'
            %% 4. Check Sign Restriction
            if all(cirf_cand(1,1:sr_horz)<=0) && all(cirf_cand(2,1:sr_horz)<=0) && all(irf_cand(3,1:sr_horz)>0) && all(irf_cand(5,1:sr_horz)<=0)
                SR_check = 1;
            elseif all(cirf_cand(1,1:sr_horz)>0) && all(cirf_cand(2,1:sr_horz)>0) && all(irf_cand(3,1:sr_horz)<0) && all(irf_cand(5,1:sr_horz)>0)
                SR_check = -1;
            else
                SR_check = 0;
            end
    end
    %% 3: If SR satisfied move on here
    if SR_check
        tel = tel+1;
        sirf(tel,:,:) = SR_check*irf_cand;
    end
end 

if count==maxcount && tel == 0
    girf =[];
    ind=0;
    return
else
    ind=1;
    girf = squeeze(nanmean(sirf));
end


