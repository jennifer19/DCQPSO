function y=f7(x)
% f1 is Schwefel function  [-500,500]
% The variabe x is a vector
%
d=length(x);
y1=sum(x.*sin(sqrt(abs(x))));
y=418.9829*d+y1;
        