function y = myfliplr(x)

%flips matrix in left-right direction i.e. changes order of the columns

y = x(:,end:-1:1);