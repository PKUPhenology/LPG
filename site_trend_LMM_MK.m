clear 
clc
data=xlsread('K:\plant phendata screen.xlsx','Sheet1','E2:N509092');
data=unique(data,'rows');
data(data(:,8)<1981,:)=[];
load 'K:\climaticfactor.mat';
P(2:end,2:end)= P(2:end,2:end)*24;
site=xlsread('K:\plant phendata screen.xlsx','Sheet1','AD2:AI457');

for si=1:length(site)
ts=mean(squeeze(T(366+60:366+151,site(si,1)+1,:)))';%spring temperature
p=sum(squeeze(P(366+60:366+151,site(si,1)+1,:)))';%spring precipitation
[b,~,~,~,stats]=regress(ts(1:35),[[1980:2014]' ones(35,1)]);
 t=mean(squeeze(T(367-61:366+59,site(si,1)+1,:)))';% temperature from last November to February
[b1,~,~,~,stats1]=regress(t(1:35),[[1980:2014]' ones(35,1)]);

    site(si,4:10)=[mean(ts) b(1) stats(3) mean(t) b1(1) stats1(3) mean(p)];
end
data(data(:,6)==4&data(:,9)<200,:)=[];
data(data(:,6)==4&data(:,9)>366,:)=[];

sptim=unique(data(:,[5 6 7]),'rows');
for i=1:length(sptim(:,1))
    [x,~]=find(data(:,5)==sptim(i,1)&data(:,6)==sptim(i,2)&data(:,7)==sptim(i,3));
    sptim(i,4)=median(data(x,9));
    sptim(i,5)=nanmin(data(x,9));
    sptim(i,6)=nanmax(data(x,9));
    data(x,11)=sptim(i,4);
    data(x,14)=sptim(i,5);
    data(x,15)=sptim(i,6);

end
data(data(:,6)==2&data(:,11)>173&(data(:,9)<1|data(:,9)>200),:)=[];%exclude abnormal values
data(data(:,10)<10,:)=[];%time series length

tims=[];
phid=[1 5; 2 7; 3 3; 4 1; 4 6];
for i=1:length(phid)
tran=unique(data(data(:,6)==phid(i,1)&data(:,7)==phid(i,2),[1 2 3 4 5 6 7 11 10 14 15]),'rows');%site-species-phen sequence
tims=[tims;tran];
end
tims(:,17)=tims(:,9);
tims(:,18)=tims(:,10);
tims(:,19)=tims(:,11);
tims(:,10:11)=nan;

tims(:,9)=tims(:,7);
tims(tims(:,6)==2&tims(:,8)<=173,9)=1;%spring flowering
tims(tims(:,6)==2&tims(:,8)>173,9)=2;%summer flowering
for i=1:size(tims,1)
     disp(['processing the ',num2str(i),'time series;'])
    trandata=data(data(:,1)==tims(i,1)&data(:,5)==tims(i,5)&data(:,6)==tims(i,6)&data(:,7)==tims(i,7),[8 9]);
    [b,~,~,~,stats]=regress(trandata(:,2),[trandata(:,1) ones(size(trandata,1),1)]);
    tims(i,10:12)=[b(1) stats(1) stats(3)];
    [sl,z,lc,uc]=mktrend(trandata(:,2),trandata(:,1));
    tims(i,13:16)=[sl z lc uc];
end


tims(tims(:,6)==4&tims(:,7)==6,7)=1;
data(:,12)=data(:,7);
data(data(:,6)==4&data(:,7)==6,12)=1;
for i=1:size(tims,1)
        data(data(:,5)==tims(i,5)&data(:,6)==tims(i,6)&data(:,7)==tims(i,7),13)=tims(i,9);
end

stim=unique(data(data(:,13)>0&data(:,10)>=10,[1:4 6 7]),'rows');
for i=1:size(stim,1)
    disp(['processing the ',num2str(i),'time series;'])
    trandata=data(data(:,1)==stim(i,1)&data(:,6)==stim(i,5)&data(:,7)==stim(i,6),[9 8 5 10]);
dat = array2table(trandata,'VariableNames',{'y','year','species','length'});
        dat.species= nominal(dat.species);
        %% Linear trend
        lme = fitlme(dat,'y~year+(1|species)+(1|length)');
        stim(i,7:9)=[nanmean(trandata(:,1)) lme.Coefficients.Estimate(2) coefTest(lme)];

        %% MK trend
        sp=unique(trandata(:,3));
        yrn=length(unique(trandata(:,2)));
        s=nan(length(sp),3+yrn*(yrn-1)/2);% Predefined a variable storing the phenological trend between each two years
    nn=0;  
        for j=1:length(sp)
                tran=trandata(trandata(:,3)==sp(j),[2 1]);
n=length(tran);
        s(j,[1 2 3])=[n  stim(j,[1 5])];
        ti=4;
        if n>=10
            nn=nn+n;
            for k=1:n-1
                for ki=k+1:n
                    if  tran(ki,1) - tran(k,1)>0&&abs(tran(ki,2) - tran(k,2))<100
                        s(j,ti) = ( tran(ki,2) - tran(k,2) ) / ( tran(ki,1) - tran(k,1) );
                    end
                    ti = ti + 1;
                end
            end
        end
    end
s(:,1:3)=[];
nc=length(s(~isnan(s)));
nn=floor((1+sqrt(1+8*nc))/2);
ss=sum(sign(s(~isnan(s))));
v = (( nn * ( nn - 1 ) * ( 2 * nn + 5 ) )) / 18;
if ss == 0
z = 0;
elseif ss > 0
z = ( ss - 1 ) / sqrt( v );
else
z = ( ss + 1 ) / sqrt( v );
end
nor = 1.96;%|z|>1.96 indicates the Sen's estimator is significant at a 95% confidence level. 
m1 = fix( ( nc - nor * sqrt( v ) ) / 2 );
m2 = fix( ( nc  + nor * sqrt( v ) ) / 2 );
s1 = sort( s (~isnan(s)));
lc = s1( m1 );% lower limit of the confidence interval
uc = s1( m2 + 1 );% upper limit of the confidence interval

sss=s1(~isnan(s1));
stim(i,10:14)=[median(sss) lc uc z length(sp)];
end
