% The funtion is used to invert the emperical CDF function f

function [xq] = ecdf_quartile(y, q)

    [f,x] = ecdf(y);
    pos = find(q<=f,1);
    xq = x(pos);