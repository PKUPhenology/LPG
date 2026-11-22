function y = photoforcing(x)
ta=x(1);
fbb=x(2);
fcc=x(3);
k=x(4);
Ft=x(5);
global t b px
bgss=b(:,2);
n=length(bgss);
pred=zeros(n,1);
for jj=1:n
    ft=0.0;
    index1=find(t(367:end,b(jj,1)-1979)>ta);
    if ~isempty(index1)
    ii=index1(1)+366;
    tn=t(:,b(jj,1)-1979);
    while ft<=Ft && ii<600
          Rf=1/(1+exp(fbb*(tn(ii)-fcc)))*(px(ii)/24).^k;
          ft=ft+Rf;
          ii=ii+1;
    end
    pred(jj)=ii-366;
    end
end
y=sqrt(mean((pred-bgss).^2)); 