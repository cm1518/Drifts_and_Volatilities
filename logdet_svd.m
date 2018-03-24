% Safe calculation of natural logarithm of determinant 
%
% function [ret] = logdet(X)
% X - matrix [n,d]
% ret - logarithm of determinant [scalar]
% 
% Written by Jason Rennie, August 2006
% Last modified: Tue Aug 22 18:11:41 2006

function [logdet] = logdet(X)
 % fn = mfilename;
 % if nargin < 1
 %   error('insufficient parameters')
 % end
  
  s = svd(X);
  logdet = sum(log(s));
  
  
  

% ChangeLog
% 8/22/06 - Use form of svd() that returns only a vector of sv's.
% 8/22/06 - First version
