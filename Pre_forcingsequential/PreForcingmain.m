clear
clc
options = optimoptions(@simulannealbnd,'Display','off');
load 'K:\climaticfactor.mat';
P(P<0.00001)=0;
P(2:end,2:end)= P(2:end,2:end)*24;
data1=importdata('K:\phendata.mat');
xlen=unique(data1(:,[1 5 6]),'rows');
phen=data1;

global t b pp
fid = fopen('K:\unichill.txt','w');
sta=[];
for i=1:size(x,1)
    disp(['Processing the ',num2str(i),' time series;']);
    index2=phen(:,1)==xlen(i,1) & phen(:,5)==xlen(i,2)& phen(:,6)==xlen(i,3);
    nn=length(nonzeros(index2));
    if nn >=10 
        t=squeeze(T(2:end,xlen(i,1)+1,:));
        pp=squeeze(P(2:end,xlen(i,1)+1,:));
        b=phen(index2,[8:9]);
        bb=b(:,2);
        mb=mean(bb);
        null=sqrt(sum((bb-mb).^2)/nn);
        fun=@preforcing;   
        x0=[-1 5 50 50];
        lb = [-100   0    0    0];
        ub = [0    max(t(:,3)) 1000 1000];
        for xunl=0:20%iteration
        [x1,fval1] = simulannealbnd(fun,x0,lb,ub,options);
        x0=x1;
        end
        nse=1-(fval1/null)^2;
        k=length(x1);
        AICc=2*nn*log(fval1)+2*nn*(k+1)/(nn-k-2);
        [co,ap,iic,pred0,predpre,predtem]=huidaipreforcing(x1);
        phen(index2,12:15)=[iic(b(:,1)-1979) pred0(b(:,1)-1979) predpre(b(:,1)-1979) predtem(b(:,1)-1979)];
        sta=[sta; [ones(39,1)*i ones(39,1)*xlen(1)  ones(39,1)*xlen(2)  ones(39,1)*xlen(3)  iic pred0 predpre predtem]];

        fprintf(fid,'%6i  %6i   %6i   %6i   %6i   %9.4f   %9.4f   %9.4f   %9.4f   %9.4f   %9.4f   %9.4f    %9.4f   %9.4f    %9.4f\r\n',i, xlen(i,1),xlen(i,2),xlen(i,3),nn,x1,null,fval1,nse,AICc,co,ap);
    else
        phen(index2,12:15)=nan;
    end
    if mod(i,5)==0|i==size(xlen,1)
        save('K:\unichillprediction.mat','phen','-mat');
        save('K:\unichillprediction8020.mat','sta','-mat');
    end
end
fclose(fid);