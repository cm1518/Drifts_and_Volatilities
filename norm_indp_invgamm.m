function [betdraw, sig2draw]=norm_indp_invgamm(y,x,mu0,v0,a,b,beta_old);
%function [betdraw, sig2draw,mustar, vstar, sig2mean, sig2var]=norminvgamm(y,x,mu0,v0,a,b);
% carry out posterior inferences of linear regression with one explanatory
% variable and indenpendent normal inverse gamma prior
n=length(y);
xx=x'*x;
xy=x'*y;
%yy=y'*y;
% conditional posterior of residual variance
astar=a+(n/2);
u=y-beta_old*x;
bstar=b+0.5*(u'*u);
tmp=gamrnd(astar,1/bstar);
sig2draw=1/tmp;
% conditional posterior of regression coefficient
vstar=1/((1/v0)+(xx/sig2draw));
mustar=vstar*((mu0/v0)+(xy/sig2draw));
betdraw=mustar+sqrt(vstar)*randn(1,1);

% sig2mean=bstar/(astar-1);
% sig2var=bstar/(((astar-1)^2)*(astar-2));



