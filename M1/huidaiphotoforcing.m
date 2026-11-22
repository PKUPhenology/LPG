function [r,p,iic,pred0,preddl,predtem] = huidaiphotoforcing(x)
ta=x(1);
fbb=x(2);
fcc=x(3);
k=x(4);
Ft=x(5);
global t b px
bgss=b(:,2);
n=length(bgss);
pred0=nan(39,1);
preddl=nan(39,1);
predtem=nan(39,1);
iic=nan(39,1);
for jj=1:39%n
    ft=0.0;
    tn=t(:,jj);
    tn80=t(:,1);
    index1=find(t(367:end,jj)>ta);
    if ~isempty(index1)
    ii=index1(1)+366;
    iic(jj,1)=ii;
    while ft<=Ft && ii<600
          Rf=1/(1+exp(fbb*(tn(ii)-fcc)))*(px(ii)/24).^k;
          ft=ft+Rf;
          ii=ii+1;
    end
    pred0(jj,1)=ii-366;
    
   
 %% daylength changes, forcing temperature no changes
iid0=iic(jj,1);
iic0=iic(1,1);
ft=0;
while ft<=Ft && iic0<600
          Rf=1/(1+exp(fbb*(tn80(iic0)-fcc)))*(px(iid0)/24).^k;
          ft=ft+Rf;
          iic0=iic0+1;
          iid0=iid0+1;
    end
    preddl(jj,1)=iic0-366;

 %% daylength no changes, forcing temperature changes
iid1=iic(1,1);
iit0=iic(jj,1);
ft=0;
while ft<=Ft && iit0<600
          Rf=1/(1+exp(fbb*(tn(iit0)-fcc)))*(px(iid1)/24).^k;
          ft=ft+Rf;
          iit0=iit0+1;
          iid1=iid1+1;
    end
    predtem(jj,1)=iit0-366;

end
end
[r,p]=corr(pred0(b(:,1)-1979),bgss);
