function [pred0,predchill,predtem] = huidaifunction(x,t)
ca=x(1);
cb=x(2);
cc=x(3);
fbb=x(4);
fcc=x(5);
Fc=x(6);
Ft=x(7);

pred0=nan(size(t,2),1);
predchill=nan(size(t,2),1);
predtem=nan(size(t,2),1);
for jj=1:size(t,2)
    ft=0.0;
    fc=0.0;
    ii=245;
    tn=t(:,jj);
    tn80=t(:,1);
    while fc<=Fc && ii<600
          Rc=1/(1+exp(ca*(tn(ii)-cc)*(tn(ii)-cc)+cb*(tn(ii)-cc)));
          fc=fc+Rc;
          ii=ii+1;
    end
    iic(jj,1)=ii;
    iit=367;
    while ft<=Ft && iit<600
          if tn(iit)<0
             Rf=0;
          else
             Rf=1/(1+exp(fbb*(tn(iit)-fcc)));
          end
          ft=ft+Rf;
          iit=iit+1;
    end
    pred0(jj,1)=max(ii,iit)-366;

   
 %% chilling changes，forcing temperature no changes
iic0=iic(jj,1);
ft=0;
iict=367;
while ft<=Ft && iict<600
          if tn80(iict)<0
             Rf=0;
          else
             Rf=1/(1+exp(fbb*(tn80(iict)-fcc)));
          end
          ft=ft+Rf;
          iict=iict+1;
    end
    predchill(jj,1)=max(iic0,iict)-366;

 %% chilling no changes，forcing temperature changes
iit0=iic(1,1);
ft=0;
iitt=367;
while ft<=Ft && iitt<600
          if tn(iitt)<0
             Rf=0;
          else
             Rf=1/(1+exp(fbb*(tn(iitt)-fcc)));
          end
          ft=ft+Rf;
          iitt=iitt+1;
    end
    predtem(jj,1)=max(iit0,iitt)-366;

end
% [r,p]=corr(pred0,bgss);
