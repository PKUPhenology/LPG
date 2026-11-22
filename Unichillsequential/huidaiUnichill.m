function [r,p,iic,pred0,predchill,predtem] = huidaiUnichill(x)
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
pred0=zeros(n,1);
predchill=zeros(n,1);
predtem=zeros(n,1);
for jj=1:n
    ft=0.0;
    fc=0.0;
    ii=245;
    tn=t(:,b(jj,1)-1979);
    tn80=t(:,1);
    while fc<=Fc && ii<600
          Rc=1/(1+exp(ca*(tn(ii)-cc)*(tn(ii)-cc)+cb*(tn(ii)-cc)));
          fc=fc+Rc;
          ii=ii+1;
    end
    iic(jj,1)=ii;
    while ft<=Ft && ii<600
          if tn(ii)<0
             Rf=0;
          else
             Rf=1/(1+exp(fbb*(tn(ii)-fcc)));
          end
          ft=ft+Rf;
          ii=ii+1;
    end
    pred0(jj,1)=ii-366;

   
 %% chilling changes, forcing temperature no changes
iic0=iic(jj,1);
ft=0;
while ft<=Ft && iic0<600
          if tn80(iic0)<0
             Rf=0;
          else
             Rf=1/(1+exp(fbb*(tn80(iic0)-fcc)));
          end
          ft=ft+Rf;
          iic0=iic0+1;
    end
    predchill(jj,1)=iic0-366;

 %% chilling no changes, forcing temperature changes
iit0=iic(1,1);
ft=0;
while ft<=Ft && iit0<600
          if tn(iit0)<0
             Rf=0;
          else
             Rf=1/(1+exp(fbb*(tn(iit0)-fcc)));
          end
          ft=ft+Rf;
          iit0=iit0+1;
    end
    predtem(jj,1)=iit0-366;

end
[r,p]=corr(pred0,bgss);
