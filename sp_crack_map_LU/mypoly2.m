function y=mypoly2(x,p)
y=nan(size(x));    
y(x<p(2))=0;
y(x>=p(2))=p(1).*(x(x>=p(2))-p(2)).^2;
end