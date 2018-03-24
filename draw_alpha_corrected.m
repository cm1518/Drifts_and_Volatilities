% Substract from the data y(t), the mean Z x B(t)
yhat = zeros(M,t);
for i = 1:t
    yhat(:,i) = y(:,i) - Z((i-1)*M+1:i*M,:)*Btdraw(:,i);
end

% This part is more tricky, check Primiceri
% Zc is a [M x M(M-1)/2] matrix defined in (A.2) page 845, Primiceri
Stemp = zeros(M,numa);
Zc = zeros(t*M,numa);
Zcblock = cell(M-1,t);
sigblockt = cell(M-1,t);
for i = 1:t
    ic = 1;
    ijc = 1;
    ikc = 1;
    for j = 2:M
        Stemp(j,((j-1)+(j-3)*(j-2)/2):ic) = - yhat(1:j-1,i)';
        ic = ic + j;
    end
    Zc((i-1)*M+1:i*M,:) =  Stemp;
    sigatemp = sigt((i-1)*M+1:i*M,:);
    siga2temp = sigatemp*sigatemp;
    for jj = 2:M
        if jj == 2
            indtemp = 1;
        else
            indtemp = 2:jj;
        end
        Zcblock{jj-1,i} = Stemp(indtemp,((jj-1)+(jj-3)*(jj-2)/2):ijc);
        ijc = ijc + jj;
    end
    for kk = 2:M
        sigblockt{kk-1,i} = siga2temp(2:kk,2:kk);
    end
end

Atdraw = [];
ind = 1;
for ii = 2:M
    Zctemp = [];
    sigmatemp = [];
    for time = 1:t
        Zctemp = [Zctemp ; Zcblock{ii-1,time}]; %#ok<AGROW>
        sigmatemp = [sigmatemp ; sigblockt{ii-1,time}]; %#ok<AGROW>
    end
    if ii == 2
        indtemp = 1;
    else
        indtemp = 2:ii;
    end
    % Draw each block of A(t)
    [Atblockdraw,log_lik2a] = carter_kohn(yhat(indtemp,:),Zctemp,sigmatemp,...
        Sblockdraw{ii-1},sizeS(ii-1),sizeS(ii-1),t,A_0_prmean(((ii-1)+(ii-3)*(ii-2)/2):ind,:),A_0_prvar(((ii-1)+(ii-3)*(ii-2)/2):ind,((ii-1)+(ii-3)*(ii-2)/2):ind));
    Atdraw = [Atdraw ; Atblockdraw]; %#ok<AGROW> % Atdraw is the final matrix of draws of A(t)
    ind = ind + ii;
end

%=====| Draw S, the covariance of A(t) (from iWishart)
% Take the SSE in the state equation of A(t)
Attemp = Atdraw(:,2:t)' - Atdraw(:,1:t-1)';
sse_2 = zeros(numa,numa);
for i = 1:t-1
    sse_2 = sse_2 + Attemp(i,:)'*Attemp(i,:);
end
% ...and subsequently draw S, the covariance matrix of A(t)
ijc = 1;
for jj=2:M
    Sinv = inv(sse_2(((jj-1)+(jj-3)*(jj-2)/2):ijc,((jj-1)+(jj-3)*(jj-2)/2):ijc) + S_prmean{jj-1});
    Sinvblockdraw = wish(Sinv,t+S_prvar(jj-1));
    Sblockdraw{jj-1} = inv(Sinvblockdraw); % this is a draw from S
    ijc = ijc + jj;
end
