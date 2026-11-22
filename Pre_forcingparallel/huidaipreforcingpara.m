function [r,p,iic,pred0,predpre,predtem] = huidaipreforcingpara(x)
fbb=x(1);
fcc=x(2);
Fc=x(3);
Ft=x(4);
global t b pp
bgss=b(:,2);
n=length(bgss);
pred0=zeros(n,1);
predpre=zeros(n,1);
predtem=zeros(n,1);
for jj=1:39
    ft=0.0;
    fc=0.0;
    ii=367;
    tn=t(:,jj);
    tn80=t(:,1);
    pn=pp(:,jj);
    while fc<=Fc && ii<600
          Rc=pn(ii);
          fc=fc+Rc;
          ii=ii+1;
    end
    iic(jj,1)=ii;
    iit=367;
    while ft<=Ft && iit<600
          Rf=1/(1+exp(fbb*(tn(iit)-fcc)));
          ft=ft+Rf;
          iit=iit+1;
    end
    pred0(jj,1)=max(ii,iit)-366;

   
 %%  precipitation changes, forcing temperature no changes
iic0=iic(jj,1);
ft=0;
iict=367;
while ft<=Ft && iict<600
          Rf=1/(1+exp(fbb*(tn80(iict)-fcc)));
          ft=ft+Rf;
          iict=iict+1;
    end
    predpre(jj,1)=max(iic0,iict)-366;

 %% precipitation no changes, forcing temperature changes
iit0=iic(1,1);
ft=0;
iitt=367;
while ft<=Ft && iitt<600
          Rf=1/(1+exp(fbb*(tn(iitt)-fcc)));
          ft=ft+Rf;
          iitt=iitt+1;
    end
    predtem(jj,1)=max(iitt,iit0)-366;

end
[r,p]=corr(pred0(b(:,1)-1979),bgss);
