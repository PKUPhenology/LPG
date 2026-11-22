clear
clc
data=xlsread('K:\plant phendata screen.xlsx','Sheet1','E2:N509092');
data(data(:,8)<1980,:)=[];
data(data(:,11)<=8,:)=[];%exclude species observed at less than 8 sites
data(data(:,6)==2&data(:,12)>173&(data(:,9)<1|data(:,9)>200),:)=[];

tims=[];
phid=[1 5; 2 7];
for i=1:size(phid,1)
    tran=unique(data(data(:,6)==phid(i,1)&data(:,7)==phid(i,2),[1 2 3 4 5 6 7 11]),'rows');%site-phen sequence
    tims=[tims;tran];
end

phen=[1 5; 2 7];%leaf unfolding, flowering
% %%  three dimension gradients based on overall data
sta=[];
stap=[];
for grd=1:size(phen,1)
    trandata=data(data(:,6)==phen(grd,1)&data(:,13)==phen(grd,2),:);
    timd=unique(trandata(:,1));
    tran1=trandata(trandata(:,8)>=1980&trandata(:,8)<=1994,:);
    ph=tran1(:,[1 2 3 4 5 12 8 9]);
    dat = array2table(ph,'VariableNames',{'siteID','lat','lon','alt','species','phenevent','yr','y'});
    dat.species = nominal(dat.species);
    dat.phenevent = nominal(dat.phenevent);
    dat.yr = nominal(dat.yr);

    lme = fitlme(dat,'y~lat+lon+alt+(1|species)+(1|yr)');%
    CI1=coefCI(lme);

    tran2=trandata(trandata(:,8)>=2000&trandata(:,8)<=2014,:);
    ph2=tran2(:,[1 2 3 4 5 12 8 9]);
    dat = array2table(ph2,'VariableNames',{'siteID','lat','lon','alt','species','phenevent','yr','y'});

    dat.species = nominal(dat.species);
    dat.phenevent = nominal(dat.phenevent);
    dat.yr = nominal(dat.yr);
    lme2 = fitlme(dat,'y~lat+lon+alt+(1|species)+(1|yr)');%
    CI2=coefCI(lme2);
    sta=[sta; phen(grd,:) lme.Coefficients.Estimate(2:4)' length(unique(ph(:,1))) lme2.Coefficients.Estimate(2:4)' length(unique(ph2(:,1))) CI1(2,:) CI2(2,:)];
    stap= [stap;phen(grd,:) lme.Coefficients.pValue(2:4)' lme2.Coefficients.pValue(2:4)'];

end

%% analysis at species level
stasps=[];
stapsps=[];
spc=[138 214 195 162 219 188 123]; %species ID
spcname={'All','PA','TM','SM','PS','UP','RP','MA'};
for i=1:length(spc)

for grd=1:2

    trandata=data(data(:,6)==phen(grd,1)&data(:,13)==phen(grd,2)&data(:,5)==spc(i),:);
    timd=unique(trandata(:,1));

tran=trandata(trandata(:,8)>=1980&trandata(:,8)<=1994,:);
tran2=trandata(trandata(:,8)>=2000&trandata(:,8)<=2014,:);
if length(unique(tran(:,1)))>=4&&length(unique(tran2(:,1)))>=4
    ph=tran(:,[1 2 3 4 5 12 8 9]);
    dat = array2table(ph,'VariableNames',{'siteID','lat','lon','alt','species','phenevent','yr','y'});

dat.species = nominal(dat.species);
dat.phenevent = nominal(dat.phenevent);
dat.yr = nominal(dat.yr);

lme = fitlme(dat,'y~lat+lon+alt+(1|yr)');%

    ph2=tran2(:,[1 2 3 4 5 12 8 9]);
    dat = array2table(ph2,'VariableNames',{'siteID','lat','lon','alt','species','phenevent','yr','y'});

dat.species = nominal(dat.species);
dat.phenevent = nominal(dat.phenevent);
dat.yr = nominal(dat.yr);

lme2 = fitlme(dat,'y~lat+lon+alt+(1|yr)');%
 CI1=coefCI(lme);
  CI2=coefCI(lme2);
stasps=[stasps; spc(i) phen(grd,:)  lme.Coefficients.Estimate(2:4)' length(unique(ph(:,1))) lme2.Coefficients.Estimate(2:4)' length(unique(ph2(:,1)))  CI1(2,:) CI2(2,:)];
stapsps= [stapsps;spc(i) phen(grd,:) lme.Coefficients.pValue(2:4)' lme2.Coefficients.pValue(2:4)'];
end
end
end
