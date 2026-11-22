function y = preforcing(x)
fbb=x(1);
fcc=x(2);
Fc=x(3);
Ft=x(4);
global t b pp
bgss=b(:,2);
n=length(bgss);
pred=zeros(n,1);
for jj=1:n
    ft=0.0;
    fc=0.0;
    ii=367;
    
    tn=t(:,b(jj,1)-1979);
    pn=pp(:,b(jj,1)-1979);
    while fc<=Fc && ii<600
          Rc=pn(ii);
          fc=fc+Rc;
          ii=ii+1;
    end
    while ft<=Ft && ii<600
          Rf=1/(1+exp(fbb*(tn(ii)-fcc)));
          ft=ft+Rf;
          ii=ii+1;
    end
    pred(jj)=ii-366;
    
end
y=sqrt(mean((pred-bgss).^2)); 