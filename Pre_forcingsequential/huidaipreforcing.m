function [r,p,iic,pred0,predpre,predtem] = huidaipreforcing(x)

fbb=x(1);
fcc=x(2);
Fc=x(3);
Ft=x(4);
global t b pp
bgss=b(:,2);
n=length(bgss);
pred0=zeros(39,1);
predpre=zeros(39,1);
predtem=zeros(39,1);
for jj=1:39%
    ft=0.0;
    fc=0.0;
    ii=367;%
    tn=t(:,jj);
    tn80=t(:,1);
    pn=pp(:,jj);
    while fc<=Fc && ii<600
          Rc=pn(ii);
          fc=fc+Rc;
          ii=ii+1;
    end
    iic(jj,1)=ii;
    while ft<=Ft && ii<600
          Rf=1/(1+exp(fbb*(tn(ii)-fcc)));
          ft=ft+Rf;
          ii=ii+1;
    end
    pred0(jj,1)=ii-366;

   
 %% precipitation changes，forcing temperature no changes
iic0=iic(jj,1);
ft=0;
while ft<=Ft && iic0<600
          Rf=1/(1+exp(fbb*(tn80(iic0)-fcc)));
          ft=ft+Rf;
          iic0=iic0+1;
    end
    predpre(jj,1)=iic0-366;

 %% precipitation no changes，forcing temperature changes
iit0=iic(1,1);
ft=0;
while ft<=Ft && iit0<600
          Rf=1/(1+exp(fbb*(tn(iit0)-fcc)));
          ft=ft+Rf;
          iit0=iit0+1;
    end
    predtem(jj,1)=iit0-366;

end
[r,p]=corr(pred0(b(:,1)-1979),bgss);
