function [bdraw,log_lik] = carter_kohn_me_with_const(y,Q,F,H,R,c,m,t,B0,V0)
% Carter and Kohn (1994), On Gibbs sampling for state space models.
% y(t) = H*b(t) + e(t)
% b(t) = c(t)+F(t)*b(t-1) + u(t)
% e(t) ~ N(0,R), where R=0 
% u(t) ~ N(0,Q(t))

% y is a Mxt data matrix, hence y(t) is Mx1


% Kalman Filter
bp = B0;
Vp = V0;
bt = zeros(t,m);
Vt = zeros(m^2,t);
log_lik = 0;
for i=1:t
    ct=c(:,i);
    Qt=Q(:,:,i);
    Ft=F(:,:,i);
    cfe = y(:,i) - H*bp;   % conditional forecast error
    K=Vp*H'/(H*Vp*H' + R); % Kalman gain
    btt=bp+K*cfe;
    Vtt=Vp-K*H*Vp;
    %f = H*Vp*H' + R;    % variance of the conditional forecast error
    %inv_f = inv(f);
    %log_lik = log_lik + log(det(f)) + cfe'*inv_f*cfe;
    %btt = bp + Vp*H'*inv_f*cfe;
    %Vtt = Vp - Vp*H'*inv_f*H*Vp;
    if i < t
        bp = ct+Ft*btt;
        Vp = Ft*Vtt*Ft' + Qt;
    end
    bt(i,:) = btt';
    Vt(:,i) = reshape(Vtt,m^2,1);
end

% draw Sdraw(T|T) ~ N(S(T|T),P(T|T))
bdraw = zeros(t,m);
bdraw(t,:) = mvnrnd(btt,Vtt,1);

% Backward recurssions
for i=1:t-1
    ct=c(:,t-i+1);
    Qt=Q(:,:,t-i+1);
    Ft=F(:,:,t-i+1);
    bf = bdraw(t-i+1,:)';
    btt = bt(t-i,:)';
    Vtt = reshape(Vt(:,t-i),m,m);
    f = Ft*Vtt*Ft' + Qt;
    inv_f = inv(f);
    cfe = bf - ct- Ft*btt;
    bmean = btt + Vtt*Ft'*inv_f*cfe;
    bvar = Vtt - Vtt*Ft'*inv_f*Ft*Vtt;
    bdraw(t-i,:) = mvnrnd(bmean,bvar,1); %bmean' + randn(1,m)*chol(bvar);
end
bdraw = bdraw';