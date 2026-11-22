function y = Unichill(x)
ca=x(1);
cb=x(2);
cc=x(3);
fbb=x(4);
fcc=x(5);
Fc=x(6);
Ft=x(7);
global t b
bgss=b(:,2);
n=length(bgss);
pred=zeros(n,1);
for jj=1:n
    ft=0.0;
    fc=0.0;
    ii=245;
    
    tn=t(:,b(jj,1)-1979);
    while fc<=Fc && ii<600
          Rc=1/(1+exp(ca*(tn(ii)-cc)*(tn(ii)-cc)+cb*(tn(ii)-cc)));
          fc=fc+Rc;
          ii=ii+1;
    end
    while ft<=Ft && ii<600
          if tn(ii)<0
             Rf=0;
          else
             Rf=1/(1+exp(fbb*(tn(ii)-fcc)));
          end
          ft=ft+Rf;
          ii=ii+1;
    end
    pred(jj)=ii-366;
    
end
y=sqrt(mean((pred-bgss).^2)); 