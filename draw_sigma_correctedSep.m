function [Sigtdraw,log_lik3,sigt,Wdraw] = draw_sigma_correctedSep(statedraw,Wdraw,yss,Zs,m_s,u2_s,M,t,sigma_prmean,sigma_prvar,W_prmean,W_prvar)
% Draw log volatilities

% In order to draw the log-volatilies, substract the mean and variance
% of the 7-component mixture of Normal approximation to the measurement
% error covariance
vart = zeros(t*M,M);
yss1 = zeros(t,M);
for i = 1:t
    for j = 1:M
        imix = statedraw(i,j);
        vart((i-1)*M+j,j) = u2_s(imix);
        yss1(i,j) = yss(i,j) - m_s(imix) + 1.2704;
    end
end


%% Old wdraw or updated???
% Sigtdraw is a draw of the diagonal log-volatilies, which will give us SIGMA(t)
[Sigtdraw,log_lik3] = carter_kohn(yss1',Zs,vart,Wdraw,M,M,t,sigma_prmean,sigma_prvar);


% Draws in Sigtdraw are in logarithmic scale (log-volatilies). Create
% original standard deviations of the VAR covariance matrix
sigtemp = eye(M);
sigt = zeros(M*t,M);
for i = 1:t
    for j = 1:M
        sigtemp(j,j) = exp(0.5*Sigtdraw(j,i));
    end
    sigt((i-1)*M+1:i*M,:) = sigtemp;
end


%=====| Draw W, the covariance of SIGMA(t) (from iWishart)
% Get first differences of Sigtdraw to compute the SSE
Sigttemp = Sigtdraw(:,2:t)' - Sigtdraw(:,1:t-1)';
sse_2 = zeros(M,M);
for i = 1:t-1
    sse_2 = sse_2 + Sigttemp(i,:)'*Sigttemp(i,:);
end
Winv = inv(sse_2 + W_prmean);
Winvdraw = wish(Winv,t+W_prvar);
Wdraw = inv(Winvdraw);  % this is a draw from W

