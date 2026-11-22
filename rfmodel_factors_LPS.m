clear
clc
data=xlsread('K:\plant phendata screen.xlsx','Sheet1','E2:N509092');
data(data(:,8)<1980,:)=[];
data(data(:,8)>2014,:)=[];
load 'K:\climaticfactor.mat';
P(P<0.00001)=0;
P(2:end,2:end)= P(2:end,2:end)*24;

sptim=unique(data(:,[5 6 7]),'rows');
for i=1:length(sptim(:,1))
    [x,~]=find(data(:,5)==sptim(i,1)&data(:,6)==sptim(i,2)&data(:,7)==sptim(i,3));
    sptim(i,4)=median(data(x,9));
    data(x,11)=length(unique(data(x,1)));
    data(x,12)=sptim(i,4);
end
data(data(:,11)<=8,:)=[];%
data(data(:,6)==2&data(:,12)>173&(data(:,9)<1|data(:,9)>200),:)=[];
data(isnan(data(:,9)),:)=[];
data(data(:,10)<10,:)=[];

%%daylength at the first day when temperature >5 °C
mask = T(368:end,:,:) > 5;
[~, first_over_T] = max(mask, [], 1);

%preseason temperature, precipitation
data8014=[];
for i=1:length(tims)
    [x,~]=find(data(:,1)==tims(i,1)&data(:,5)==tims(i,5)&data(:,6)==tims(i,6)&data(:,7)==tims(i,7));
    phtims=data(x,[8 9]);
    mph=ceil(mean(phtims(:,2)));
    data(x,14)=mph;
    data(x,15)=phtims(:,2)-mph;


    %% Forcing temperature related variables
    for mlp=1:18
        t0=mean(squeeze(T(367+mph-2-mlp*5:367+mph+2,tims(i,1)+1,:)))';
        [r(mlp),~]=corr(phtims(:,2),t0(phtims(:,1)-1979));
    end
    [~,loc]=max(abs(r));%
    tp=nanmean(squeeze(T(367+mph-2-loc*5:367+mph+2,tims(i,1)+1,:)))';%preseason temperature
    ts_sd=nanstd(squeeze(T(367+60:367+150,tims(i,1)+1,:))-nanmean(squeeze(T(367+60:367+150,tims(i,1)+1,:)),2))';%sd for spring temperature
Rt(i)=r(loc);
Loct(i)=loc*5+5;
    %% chilling exposure  related variables
    t=squeeze(T(368-61-61:367+mph,tims(i,1)+1,:));
    cdays=nansum(t<=8&t>=0)';%
    tw=nanmean(squeeze(T(368-31:367+60,tims(i,1)+1,:)))';%winter temperature
    tw_sd=nanstd(squeeze(T(368-31:367+60,tims(i,1)+1,:))-nanmean(squeeze(T(368-31:367+60,tims(i,1)+1,:)),2))';%sd for winter temperature


    %% precipitation  related variables

    for mlp=1:18
        p=sum(squeeze(P(367+mph-2-mlp*5:367+mph+2,tims(i,1)+1,:)))';
        [r(mlp),~]=corr(phtims(:,2),p(phtims(:,1)-1979));
    end
    [~,loc]=max(abs(r));%
    pp=nansum(squeeze(P(367+mph-2-loc*5:367+mph+2,tims(i,1)+1,:)))';%preseason precipitation
    pps=nansum(squeeze(P(367+60:367+150,tims(i,1)+1,:)))';%spring precipitation
Rp(i)=r(loc);
Locp(i)=loc*5+5;

    %%  daylength related variables
    lat=tims(i,2);%latitude
    doy=1:365;%
    ao=asin(0.39795*cos(0.2163108+2*atan(0.9671396*tan(0.00860*(doy-186)))));
    px=24-24*acos((sin(0.833*pi/180)+sin(lat*pi/180)*sin(ao))./(cos(lat*pi/180)*cos(ao)))/pi;%
    DL=px(squeeze(first_over_T(1,tims(i,1)+1,:)))';
    DLm=px(mph)*ones(39,1);
    DLm0=px(phtims(1,2));

    %% background variables
    tpa=nanmean(squeeze(T(368:733,tims(i,1)+1,:)))';%annual temperature
    tpa_range=nanmax(squeeze(T(368:733,tims(i,1)+1,:)))'-nanmin(squeeze(T(368:733,tims(i,1)+1,:)))';%annual temperature range
    ppa=nansum(squeeze(P(368:733,tims(i,1)+1,:)))';%annual precipitation
    px_range(1:39,1)=nanmax(px)-nanmin(px);%annual range in daylength
    tran8014=[tp-mean(tp)  ts_sd cdays-mean(cdays) tw-mean(tw) tw_sd (pp-min(pp))/(max(pp)-min(pp))   DL-mean(DL) DLm tpa-mean(tpa) tpa_range (ppa-mean(ppa))/(max(ppa)-min(ppa)) px_range];%tcd-mean(tcd)%ts-mean(ts), (pps-min(pps))/(max(pps)-min(pps))
    tran8014=[tran8014 [tran8014(1,1:10) DLm0].*ones(39,1)];
    data(x,16:38)=tran8014(phtims(:,1)-1979,:);%16-18 forcing, 19-21 chilling, 22-23 prec, 24-25 photoperiod, 26-29 background
    data8014=[data8014; [tims(i,1:7) data(x(1),13:14)].*ones(39,1) (1980:2018)' tran8014];%11-13 forcing, 14-16 chilling, 17-18 prec, 19-20 photoperiod, 21-24 background
    
end


%% random forest model
for ph=1:2
    for regionID=1:4
        tData=data(data(:,6)==ph&data(:,13)==regionID,[16:27]); 
        speciesID=categorical(data(data(:,6)==ph&data(:,13)==regionID,5));
        envirData=array2table(tData,'VariableNames',{'tp', 'ts_sd', 'cdays', 'tw', 'tw_sd', 'pp',  'DL', 'DLm', 'tpa', 'tpa_range', 'ppa', 'px_range'});%'tcd', ,'lat','lon','alt' 'ts','pps',
        spData=table(speciesID,'VariableNames',{'species'});
        xData=[envirData spData];
        yData=data(data(:,6)==ph&data(:,13)==regionID,15);
        cv=cvpartition(size(xData,1),'HoldOut',0.3);
        idx=cv.test;
        xTrain=xData(~idx,:);
        yTrain=yData(~idx,:);
        xTest=xData(idx,:);
        yTest=yData(idx,:);
        
        numTrees=500;
        rng(42);
        rfmodel=TreeBagger(numTrees,xTrain,yTrain,'Method','regression','OOBPrediction','on','OOBPredictorImportance','on','CategoricalPredictors','species');
        eval(['rfmodel',num2str((ph-1)*4+regionID),'=rfmodel;']);
        yPred=predict(rfmodel,xTest);
        data(data(:,6)==ph&data(:,13)==regionID,41)=predict(rfmodel,xData)+data(data(:,6)==ph&data(:,13)==regionID,14);%
        %model performance
        rmse(regionID,ph)=sqrt(mean((yTest-yPred).^2));
        rR=corrcoef(yTest,yPred);
        r2(regionID,ph)=rR(1,2).^2;
        yPredtrain=predict(rfmodel,xTrain);
        rmsetrain(regionID,ph)=sqrt(mean((yTrain-yPredtrain).^2));
        rRt=corrcoef(yTrain,yPredtrain);
        r2train(regionID,ph)=rRt(1,2).^2;
        importance(regionID+(ph-1)*4,:)=rfmodel.OOBPermutedPredictorDeltaError;
        disp(['feature importance：',num2str(importance(regionID,:))]);
        tData8014=data8014(data8014(:,6)==ph&data8014(:,8)==regionID,[11:22]); 
        speciesID8014=categorical(data8014(data8014(:,6)==ph&data8014(:,8)==regionID,5));
        envirData=array2table(tData8014,'VariableNames',{'tp', 'ts_sd', 'cdays',  'tw', 'tw_sd', 'pp', 'DL', 'DLm', 'tpa', 'tpa_range', 'ppa', 'px_range'});%'tcd',,'lat','lon','alt' 'ts', 'pps',
        spData=table(speciesID8014,'VariableNames',{'species'});
        xData8014=[envirData spData];
        data8014(data8014(:,6)==ph&data8014(:,8)==regionID,36)=predict(rfmodel,xData8014)+data8014(data8014(:,6)==ph&data8014(:,8)==regionID,9);%

        for f=1:4
            cfactor=[1 0 2 0 0 3 4 0 0 0 0 0;31:38 27:30];
            cFactor=[16:27];
            cFactor(cfactor(1,:)==f)=cfactor(2,cfactor(1,:)==f);
            tData1=data(data(:,6)==ph&data(:,13)==regionID,[cFactor]);
            envirData=array2table(tData1,'VariableNames',{'tp',  'ts_sd', 'cdays',  'tw', 'tw_sd', 'pp', 'DL', 'DLm', 'tpa', 'tpa_range', 'ppa', 'px_range'});
            spData=table(speciesID,'VariableNames',{'species'});
            xTest1=[envirData spData];
            yTest1=data(data(:,6)==ph&data(:,13)==regionID,15);
            yPred1=predict(rfmodel,xTest1);
            rmse1(regionID,f)=sqrt(mean((yTest1-yPred1).^2));
            rtran=corrcoef(yTest1,yPred1);
            r21(regionID,f)=rtran(1,2).^2;
            data(data(:,6)==ph&data(:,13)==regionID,41+f)=yPred1+data(data(:,6)==ph&data(:,13)==regionID,14);%

            tData18014=data8014(data8014(:,6)==ph&data8014(:,8)==regionID,[cFactor-5]);
            speciesID18014=categorical(data8014(data8014(:,6)==ph&data8014(:,8)==regionID,5));
            envirData=array2table(tData18014,'VariableNames',{'tp',  'ts_sd', 'cdays', 'tw', 'tw_sd', 'pp',  'DL', 'DLm', 'tpa', 'tpa_range', 'ppa', 'px_range'}); 
            spData=table(speciesID18014,'VariableNames',{'species'});
            xTest18014=[envirData spData];
            yPred2=predict(rfmodel,xTest18014);
            data8014(data8014(:,6)==ph&data8014(:,8)==regionID,36+f)=yPred2+data8014(data8014(:,6)==ph&data8014(:,8)==regionID,9);
        end
    end
    rfTims=unique(data(data(:,6)==ph,1:7),'rows');
    rfTims1=rfTims;
    rfTims2=rfTims;
    for i=1:size(rfTims,1)
        phTims=[];
        phTims1=[];
        phTims=data(data(:,1)==rfTims(i,1)&data(:,5)==rfTims(i,5)&data(:,6)==ph,[8 42:45 41 16]);
        for j=1:5
            if ~isempty(phTims)&&size(phTims,1)>=5
                [b,bint,r,rint,stats]=regress(phTims(:,j+1),[phTims(:,7) ones(size(phTims,1),1)]);
                rfTims(i,8+(j-1)*3:10+(j-1)*3)=[b(1) stats(1) stats(3)];
            else
                rfTims(i,8+(j-1)*3:10+(j-1)*3)=nan;
            end
        end
        [b,bint,r,rint,stats]=regress(phTims(:,6),[phTims(:,1) ones(size(phTims,1),1)]);
        rfTims2(i,8:10)=[b(1) stats(1) stats(3)];
    end
    rfTims(isnan(rfTims(:,22)),:)=[];
    if ph==1
        rfTims0=rfTims;
    end
    staslope=[];
    %% 计算modeled LPS
    year=[1980 2000];
    trandata=data8014(data8014(:,6)==ph,:);
    for yr=1:2
        for j=1:5
            phlps=data8014(data8014(:,6)==ph&data8014(:,10)>=year(yr)&data8014(:,10)<=year(yr)+14,[1 2 3 4 5 9 10 35+j]);
            datlps = array2table(phlps,'VariableNames',{'siteID','lat','lon','alt','species','phenevent','yyr','y'});
            datlps.species = nominal(datlps.species);
            datlps.phenevent = nominal(datlps.phenevent);
            datlps.yr = nominal(datlps.yyr);
            lme = fitlme(datlps,'y~lat+alt+lon+(1|species)+(1|yyr)');
            staslope=[staslope; ph lme.Coefficients.Estimate(2)' lme.Coefficients.pValue(2)'];
        end
    end
    if ph==1
        staslope0=staslope;
    end

    