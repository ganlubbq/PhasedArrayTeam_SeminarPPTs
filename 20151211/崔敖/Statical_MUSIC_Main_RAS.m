clear;
for STCounter=1:2
    OriginalRightRate=0;
    RACRightRate=0;
    KSRightRate=0;
    MSRightRate=0;
    KKRightRate=0;
%clc;
%% Init and Settings
RotateAngle=83:0.1:90;
for AngleCounter=1:1:length(RotateAngle)
  OrignalRight=0;
  RACRight=0;
  KSRight=0;
  MSRight=0;
  KKRight=0;
  StaticCounterNum=200;
  tic;
  for StaticCounter=1:StaticCounterNum
    bb=double.empty();
    RMatrixTime=double.empty();
    ExpSetting=struct('SampleFrequency',48000, ...
                        'WaveSpeed',340, ...
                        'SNR',10, ...
                        'TimeSnaps',200);

%     SignalSetting=struct( 'SignalNum',2, ...
%                              'SignalFrequency',[1000,1000], ...
%                              'IncidentAngle',[round(unifrnd(10,80)),round(unifrnd(100,170))], ...
%                              'Power',[0,0]);
                         
%     SignalSetting=struct( 'SignalNum',1, ...
%                              'SignalFrequency',[1000], ...
%                              'IncidentAngle',[50,], ...
%                              'Power',[0]);
                         
    SignalSetting=struct( 'SignalNum',2, ...
                             'SignalFrequency',[1000,1000], ...
                             'IncidentAngle',[34.5,48], ...
                             'Power',[0,0]);

    ArraySetting=struct(  'ElementNum',5, ...
                            'ElementSpacingTime',2, ...
                            'AssistArrayElementSpacingTime',2,...
                            'AssistArrayRotateAngle',RotateAngle(AngleCounter));

    Signal=struct(        'SignalNum',SignalSetting.SignalNum, ...
                            'DigitalFrequency',SignalSetting.SignalFrequency.*2*pi./ExpSetting.SampleFrequency, ...
                            'WaveLength',ExpSetting.WaveSpeed./SignalSetting.SignalFrequency, ...
                            'IncidentAngle',SignalSetting.IncidentAngle, ...
                            'RotateIncidentAngle_P',double.empty, ...
                            'RotateIncidentAngle_N',double.empty, ...
                            'SteeringSpacing',double.empty(), ...
                            'SteeringVectors',double.empty(), ...
                            'AssistSteeringSpacing_P',double.empty(), ...
                            'AssistSteeringVectors_P',double.empty(), ...
                            'AssistSteeringSpacing_N',double.empty(), ...
                            'AssistSteeringVectors_N',double.empty());
%(ArraySetting.ElementNum-1)*ArraySetting.ElementSpacingTime/0.5+1
    Array=struct(         'ElementNum',(ArraySetting.ElementNum-1)*ArraySetting.ElementSpacingTime/0.5+1, ...
                            'ElementSpacing',min(Signal.WaveLength)/2, ...
                            'SteeringAngles',0:0.5:180, ...
                            'SteeringSpacing',double.empty(), ...
                            'SteeringVectors',double.empty(), ...
                            'MusicSpectrum',double.empty(), ...
                            'RxSignal', double.empty(), ...
                            'RMatrixTime',double.empty(), ...
                            'RMatrix',double.empty(), ...
                            'NoiseBase',double.empty(), ...
                            'MusicSpectrumMean',double.empty(), ...
                            'MusicVar',double.empty(), ...
                            'GroupDelay',double.empty(), ...
                            'MGDSpectrum',double.empty(), ...
                            'MusicBound',double.empty());

    AssistArray_P=struct(   'ElementNum',ArraySetting.ElementNum, ...
                            'ElementSpacing',ArraySetting.AssistArrayElementSpacingTime*min(Signal.WaveLength), ...
                            'RotateAngle',ArraySetting.AssistArrayRotateAngle, ...
                            'SteeringAngles',0:0.5:180, ...
                            'SteeringSpacing',double.empty(), ...
                            'SteeringVectors',double.empty(), ...
                            'MusicSpectrum',double.empty(), ...
                            'RxSignal', double.empty(), ...
                            'RMatrixTime',double.empty(), ...
                            'RMatrix',double.empty(), ...
                            'NoiseBase',double.empty(), ...
                            'GroupDelay',double.empty(), ...
                            'MGDSpectrum',double.empty());

    AssistArray_N=struct(   'ElementNum',ArraySetting.ElementNum, ...
                            'ElementSpacing',ArraySetting.AssistArrayElementSpacingTime*min(Signal.WaveLength), ...
                            'RotateAngle',ArraySetting.AssistArrayRotateAngle, ...
                            'SteeringAngles',0:0.5:180, ...
                            'SteeringSpacing',double.empty(), ...
                            'SteeringVectors',double.empty(), ...
                            'MusicSpectrum',double.empty(), ...
                            'RxSignal', double.empty(), ...
                            'RMatrixTime',double.empty(), ...
                            'RMatrix',double.empty(), ...
                            'NoiseBase',double.empty(), ...
                            'GroupDelay',double.empty(), ...
                            'MGDSpectrum',double.empty());
    tic
      Signal.RotateIncidentAngle_P=SignalSetting.IncidentAngle+AssistArray_P.RotateAngle;
      Signal.RotateIncidentAngle_N=SignalSetting.IncidentAngle-AssistArray_N.RotateAngle;
      Array.SteeringSpacing=2*pi*SignalSetting.SignalFrequency(1).*Array.ElementSpacing.*cosd(Array.SteeringAngles)./ExpSetting.WaveSpeed;
      for Counter=1:1:Array.ElementNum
        Array.SteeringVectors(Counter,:)=exp(-1i*(Counter-1).*Array.SteeringSpacing);
      end

      AssistArray_P.SteeringSpacing=2*pi*SignalSetting.SignalFrequency(1).*AssistArray_P.ElementSpacing.*cosd(AssistArray_P.SteeringAngles+ArraySetting.AssistArrayRotateAngle)./ExpSetting.WaveSpeed;
      for Counter=1:1:AssistArray_P.ElementNum
        AssistArray_P.SteeringVectors(Counter,:)=exp(-1i*(Counter-1).*AssistArray_P.SteeringSpacing);
      end

      AssistArray_N.SteeringSpacing=2*pi*SignalSetting.SignalFrequency(1).*AssistArray_N.ElementSpacing.*cosd(AssistArray_N.SteeringAngles-ArraySetting.AssistArrayRotateAngle)./ExpSetting.WaveSpeed;
      for Counter=1:1:AssistArray_N.ElementNum
        AssistArray_N.SteeringVectors(Counter,:)=exp(-1i*(Counter-1).*AssistArray_N.SteeringSpacing);
      end

      Signal.SteeringSpacing=2*pi*SignalSetting.SignalFrequency(1).*Array.ElementSpacing.*cosd(Signal.IncidentAngle)./ExpSetting.WaveSpeed;
      for Counter=1:1:Array.ElementNum
        Signal.SteeringVectors(Counter,:)=exp(-1i*(Counter-1).*Signal.SteeringSpacing);
      end

      Signal.AssistSteeringSpacing_P=2*pi*SignalSetting.SignalFrequency(1).*AssistArray_P.ElementSpacing.*cosd(Signal.RotateIncidentAngle_P)./ExpSetting.WaveSpeed;
      for Counter=1:1:AssistArray_P.ElementNum
        Signal.AssistSteeringVectors_P(Counter,:)=exp(-1i*(Counter-1).*Signal.AssistSteeringSpacing_P);
      end

      Signal.AssistSteeringSpacing_N=2*pi*SignalSetting.SignalFrequency(1).*AssistArray_N.ElementSpacing.*cosd(Signal.RotateIncidentAngle_N)./ExpSetting.WaveSpeed;
      for Counter=1:1:AssistArray_N.ElementNum
        Signal.AssistSteeringVectors_N(Counter,:)=exp(-1i*(Counter-1).*Signal.AssistSteeringSpacing_N);
      end

      for Counter=1:1:ExpSetting.TimeSnaps
        Array.RxSignal(:,Counter)=Signal.SteeringVectors*sqrt(2)*[cos(Signal.DigitalFrequency(1)*Counter);sin(Signal.DigitalFrequency(1)*Counter)];
      end

      for Counter=1:1:Array.ElementNum
        Array.RxSignal(Counter,:)=awgn(Array.RxSignal(Counter,:),ExpSetting.SNR,'measured');
      end

      for Counter=1:1:ExpSetting.TimeSnaps
        AssistArray_P.RxSignal(:,Counter)=Signal.AssistSteeringVectors_P*sqrt(2)*[cos(Signal.DigitalFrequency(1)*Counter);sin(Signal.DigitalFrequency(1)*Counter)];
      end
      for Counter=1:1:AssistArray_P.ElementNum
        AssistArray_P.RxSignal(Counter,:)=awgn(AssistArray_P.RxSignal(Counter,:),ExpSetting.SNR,'measured');
      end

      for Counter=1:1:ExpSetting.TimeSnaps
        AssistArray_N.RxSignal(:,Counter)=Signal.AssistSteeringVectors_N*sqrt(2)*[cos(Signal.DigitalFrequency(1)*Counter);sin(Signal.DigitalFrequency(1)*Counter)];
      end
      for Counter=1:1:AssistArray_N.ElementNum
        AssistArray_N.RxSignal(Counter,:)=awgn(AssistArray_N.RxSignal(Counter,:),ExpSetting.SNR,'measured');
      end

      for Counter=1:1:ExpSetting.TimeSnaps
        Array.RMatrixTime(:,:,Counter)=Array.RxSignal(:,Counter)*(Array.RxSignal(:,Counter))';
      end
      Array.RMatrix=sum(Array.RMatrixTime,3)/ExpSetting.TimeSnaps;

      for Counter=1:1:ExpSetting.TimeSnaps
        AssistArray_P.RMatrixTime(:,:,Counter)=AssistArray_P.RxSignal(:,Counter)*(AssistArray_P.RxSignal(:,Counter))';
      end
      AssistArray_P.RMatrix=sum(AssistArray_P.RMatrixTime,3)/ExpSetting.TimeSnaps;

      for Counter=1:1:ExpSetting.TimeSnaps
        AssistArray_N.RMatrixTime(:,:,Counter)=AssistArray_N.RxSignal(:,Counter)*(AssistArray_N.RxSignal(:,Counter))';
      end
      AssistArray_N.RMatrix=sum(AssistArray_N.RMatrixTime,3)/ExpSetting.TimeSnaps;


      [U,S,V]=svd(Array.RMatrix);
      Array.NoiseBase=U(:,Signal.SignalNum+1:end);

      [U,S,V]=svd(AssistArray_P.RMatrix);
      AssistArray_P.NoiseBase=U(:,Signal.SignalNum+1:end);

      [U,S,V]=svd(AssistArray_N.RMatrix);
      AssistArray_N.NoiseBase=U(:,Signal.SignalNum+1:end);

      for Counter=1:1:size(Array.SteeringAngles,2)
        Array.MusicSpectrum(Counter)=abs(inv((Array.SteeringVectors(:,Counter))'*Array.NoiseBase*Array.NoiseBase'*Array.SteeringVectors(:,Counter)));
      end
      Array.MusicSpectrum=Array.MusicSpectrum./max(Array.MusicSpectrum);
      
      Array.GroupDelay=0;
      for Counter=1:size(Array.NoiseBase,2)
          Array.GroupDelay=Array.GroupDelay+angle((Array.NoiseBase(:,Counter))'*Array.SteeringVectors);
      end
      Array.MGDSpectrum=[0,diff(Array.GroupDelay)].*Array.MusicSpectrum;
      Array.MGDSpectrum=abs(Array.MGDSpectrum./max(Array.MGDSpectrum));

      for Counter=1:1:size(AssistArray_P.SteeringAngles,2)
        AssistArray_P.MusicSpectrum(Counter)=abs(inv((AssistArray_P.SteeringVectors(:,Counter))'*AssistArray_P.NoiseBase*AssistArray_P.NoiseBase'*AssistArray_P.SteeringVectors(:,Counter)));
      end
      
      AssistArray_P.GroupDelay=0;
      for Counter=1:size(AssistArray_P.NoiseBase,2)
          AssistArray_P.GroupDelay=AssistArray_P.GroupDelay+angle((AssistArray_P.NoiseBase(:,Counter))'*AssistArray_P.SteeringVectors);
      end
      AssistArray_P.MGDSpectrum=[0,diff(AssistArray_P.GroupDelay)].*AssistArray_P.MusicSpectrum;
      AssistArray_P.MGDSpectrum=abs(AssistArray_P.MGDSpectrum./max(AssistArray_P.MGDSpectrum));

      for Counter=1:1:size(AssistArray_N.SteeringAngles,2)
        AssistArray_N.MusicSpectrum(Counter)=abs(inv((AssistArray_N.SteeringVectors(:,Counter))'*AssistArray_N.NoiseBase*AssistArray_N.NoiseBase'*AssistArray_N.SteeringVectors(:,Counter)));
      end
      
      AssistArray_N.GroupDelay=0;
      for Counter=1:size(AssistArray_N.NoiseBase,2)
          AssistArray_N.GroupDelay=AssistArray_N.GroupDelay+angle((AssistArray_N.NoiseBase(:,Counter))'*AssistArray_N.SteeringVectors);
      end
      AssistArray_N.MGDSpectrum=[0,diff(AssistArray_N.GroupDelay)].*AssistArray_N.MusicSpectrum;
      AssistArray_N.MGDSpectrum=abs(AssistArray_N.MGDSpectrum./max(AssistArray_N.MGDSpectrum));
      
      
      ResultSpectrum=(AssistArray_P.MusicSpectrum.*AssistArray_N.MusicSpectrum).^(1/2);%.*Array.MusicSpectrum).^(1/3);
      ResultSpectrum=ResultSpectrum./max(ResultSpectrum);
      
      ResultMGD=(AssistArray_P.MGDSpectrum.*AssistArray_N.MGDSpectrum);
      ResultMGD=ResultMGD./max(ResultMGD);
      
      kk=(AssistArray_P.MusicSpectrum.*AssistArray_N.MusicSpectrum.*ResultMGD).^(1/3);
      kk=kk./max(kk);
      
      Spacing1=2*pi*SignalSetting.SignalFrequency(1).*AssistArray_P.ElementSpacing.*cosd(ArraySetting.AssistArrayRotateAngle)*cosd(AssistArray_P.SteeringAngles)./ExpSetting.WaveSpeed;
      Spacing2=2*pi*SignalSetting.SignalFrequency(1).*AssistArray_P.ElementSpacing.*sind(ArraySetting.AssistArrayRotateAngle)*sind(AssistArray_P.SteeringAngles)./ExpSetting.WaveSpeed;
for Ang=1:1:size(Spacing1,2)
    for Counter=1:1:AssistArray_P.ElementNum
        bb(Counter,Ang)=2*cos((Counter-1).*Spacing2(:,Ang))*exp(-1i*(Counter-1).*Spacing1(:,Ang));
    end
end

Rx=AssistArray_P.RxSignal+AssistArray_N.RxSignal;
for Counter=1:1:ExpSetting.TimeSnaps
    RMatrixTime(:,:,Counter)=Rx(:,Counter)*(Rx(:,Counter))';
end
RMtarix=sum(RMatrixTime,3)/ExpSetting.TimeSnaps;
[q,w,e]=svd(RMtarix);
NoiseBase=q(:,Signal.SignalNum+1:end);
for Counter=1:1:size(bb,2)
    ms(Counter)=abs(inv((bb(:,Counter))'*NoiseBase*NoiseBase'*bb(:,Counter)));
end
%ms=ms./max(ms);
ks=(AssistArray_P.MusicSpectrum.*AssistArray_N.MusicSpectrum.*ms).^(1/3);
ks=ks./max(ks);

% figure(1);
% subplot(2,1,1);
% plot(0:0.5:180,10*log10(Array.MusicSpectrum),'b',0:0.5:180,10*log10(ks),'r');
% set(gca,'XLim',[0,180],'XTick',0:10:180,'YLim',[-60,0],'YTick',-60:10:0);
% 
% subplot(2,1,2);
% plot(0:0.5:180,10*log10(ResultSpectrum),'b',0:0.5:180,10*log10(ms),'g',0:0.5:180,10*log10(ks),'r');
% set(gca,'XLim',[0,180],'XTick',0:10:180,'YLim',[-60,0],'YTick',-60:10:0);
      
      
      
%       figure(1)
%       plot(Array.SteeringAngles,Array.MusicSpectrum,'r:',Array.SteeringAngles,ResultSpectrum,'b');
%       set(gca,'Xlim',[0,180],'XTick',0:5:180);
%       figure(1)
%       title('Result Spectrums');
%       plot(0:0.5:180,kk,'r');
%       set(gca,'XTick',0:10:180);
%       set(gca,'XMinorTick','on');
%             
%       figure(2)
%       title('ORG Spectrums');
%       plot(0:0.5:180,[0,diff(AssistArray_P.GroupDelay)].*AssistArray_P.MusicSpectrum,'b',0:0.5:180,[0,diff(AssistArray_N.GroupDelay)].*AssistArray_N.MusicSpectrum,'r');
%       set(gca,'XTick',0:10:180);
%       set(gca,'XMinorTick','on');
%       
%       figure(3)
%       title('MUL MGD');
%       plot(0:0.5:180,[0,diff(AssistArray_P.GroupDelay)].*[0,diff(AssistArray_N.GroupDelay)]);
%       set(gca,'XTick',0:10:180);
%       set(gca,'XMinorTick','on');
      
      [CheckOrg,CheckOrgIndex]=findpeaks(Array.MusicSpectrum);
      [~,Peak1]=max(CheckOrg);
      CheckOrg(Peak1)=CheckOrg(Peak1)-CheckOrg(Peak1);
      [~,Peak2]=max(CheckOrg);
      if (abs(Array.SteeringAngles(CheckOrgIndex(Peak1))-SignalSetting.IncidentAngle(1))<=5 || abs(Array.SteeringAngles(CheckOrgIndex(Peak1))-SignalSetting.IncidentAngle(2))<=5) && (abs(Array.SteeringAngles(CheckOrgIndex(Peak2))-SignalSetting.IncidentAngle(1))<=5 || abs(Array.SteeringAngles(CheckOrgIndex(Peak2))-SignalSetting.IncidentAngle(2))<=5)
          OrignalRight=OrignalRight+1;
      end
      
      [CheckImp,CheckImpIndex]=findpeaks(ResultSpectrum);
      [~,Peak1]=max(CheckImp);
      CheckImp(Peak1)=CheckImp(Peak1)-CheckImp(Peak1);
      [~,Peak2]=max(CheckImp);
      if (abs(Array.SteeringAngles(CheckImpIndex(Peak1))-SignalSetting.IncidentAngle(1))<=5 || abs(Array.SteeringAngles(CheckImpIndex(Peak1))-SignalSetting.IncidentAngle(2))<=5) && (abs(Array.SteeringAngles(CheckImpIndex(Peak2))-SignalSetting.IncidentAngle(1))<=5 || abs(Array.SteeringAngles(CheckImpIndex(Peak2))-SignalSetting.IncidentAngle(2))<=5)
          RACRight=RACRight+1;
      end
      
      [CheckImp,CheckImpIndex]=findpeaks(ks);
      [~,Peak1]=max(CheckImp);
      CheckImp(Peak1)=CheckImp(Peak1)-CheckImp(Peak1);
      [~,Peak2]=max(CheckImp);
      if (abs(Array.SteeringAngles(CheckImpIndex(Peak1))-SignalSetting.IncidentAngle(1))<=5 || abs(Array.SteeringAngles(CheckImpIndex(Peak1))-SignalSetting.IncidentAngle(2))<=5) && (abs(Array.SteeringAngles(CheckImpIndex(Peak2))-SignalSetting.IncidentAngle(1))<=5 || abs(Array.SteeringAngles(CheckImpIndex(Peak2))-SignalSetting.IncidentAngle(2))<=5)
          KSRight=KSRight+1;
      end
      
      [CheckImp,CheckImpIndex]=findpeaks(ms);
      [~,Peak1]=max(CheckImp);
      CheckImp(Peak1)=CheckImp(Peak1)-CheckImp(Peak1);
      [~,Peak2]=max(CheckImp);
      if (abs(Array.SteeringAngles(CheckImpIndex(Peak1))-SignalSetting.IncidentAngle(1))<=5 || abs(Array.SteeringAngles(CheckImpIndex(Peak1))-SignalSetting.IncidentAngle(2))<=5) && (abs(Array.SteeringAngles(CheckImpIndex(Peak2))-SignalSetting.IncidentAngle(1))<=5 || abs(Array.SteeringAngles(CheckImpIndex(Peak2))-SignalSetting.IncidentAngle(2))<=5)
          MSRight=MSRight+1;
      end
      
      [CheckImp,CheckImpIndex]=findpeaks(kk);
      [~,Peak1]=max(CheckImp);
      CheckImp(Peak1)=CheckImp(Peak1)-CheckImp(Peak1);
      [~,Peak2]=max(CheckImp);
      if (abs(Array.SteeringAngles(CheckImpIndex(Peak1))-SignalSetting.IncidentAngle(1))<=5 || abs(Array.SteeringAngles(CheckImpIndex(Peak1))-SignalSetting.IncidentAngle(2))<=5) && (abs(Array.SteeringAngles(CheckImpIndex(Peak2))-SignalSetting.IncidentAngle(1))<=5 || abs(Array.SteeringAngles(CheckImpIndex(Peak2))-SignalSetting.IncidentAngle(2))<=5)
          KKRight=KKRight+1;
      end
%       CheckOrg=Array.MusicSpectrum;
%       [~,Peak1]=max(CheckOrg);
%       CheckOrg(Peak1)=CheckOrg(Peak1)-CheckOrg(Peak1);
%       [~,Peak2]=max(CheckOrg);
%       if (Array.SteeringAngles(Peak1)==SignalSetting.IncidentAngle(1) || Array.SteeringAngles(Peak1)==SignalSetting.IncidentAngle(2)) && (Array.SteeringAngles(Peak2)==SignalSetting.IncidentAngle(1) || Array.SteeringAngles(Peak2)==SignalSetting.IncidentAngle(2))
%           OrignalRight=OrignalRight+1;
%       end
%       
%       CheckImp=ResultSpectrum;
%       [~,Peak1]=max(CheckImp);
%       CheckImp(Peak1)=CheckImp(Peak1)-CheckImp(Peak1);
%       [~,Peak2]=max(CheckImp);
%       if (Array.SteeringAngles(Peak1)==SignalSetting.IncidentAngle(1) || Array.SteeringAngles(Peak1)==SignalSetting.IncidentAngle(2)) && (Array.SteeringAngles(Peak2)==SignalSetting.IncidentAngle(1) || Array.SteeringAngles(Peak2)==SignalSetting.IncidentAngle(2))
%           RACRight=RACRight+1;
%       end
  end
  toc;
  %Res(AngleCounter,:)=ResultSpectrum;
  disp(RotateAngle(AngleCounter));
  OriginalRightRate(AngleCounter)=OrignalRight/StaticCounterNum;
  RACRightRate(AngleCounter)=RACRight/StaticCounterNum;
  KSRightRate(AngleCounter)=KSRight/StaticCounterNum;
  MSRightRate(AngleCounter)=MSRight/StaticCounterNum;
  KKRightRate(AngleCounter)=KKRight/StaticCounterNum;
  disp(OriginalRightRate(AngleCounter));
  disp(RACRightRate(AngleCounter));
  disp(KSRightRate(AngleCounter));
  disp(MSRightRate(AngleCounter));
  disp(KKRightRate(AngleCounter));
end
OriginalRightRateTotal(STCounter,:)=OriginalRightRate;
RACRightRateTotal(STCounter,:)=RACRightRate;
KSRightRateTotal(STCounter,:)=KSRightRate;
MSRightRateTotal(STCounter,:)=MSRightRate;
KKRightRateTotal(STCounter,:)=KKRightRate;

OriginalRightRateMean=mean(OriginalRightRateTotal);
RACRightRateMean=mean(RACRightRateTotal);
KSRightRateMean=mean(KSRightRateTotal);
MSRightRateMean=mean(MSRightRateTotal);
KKRightRateMean=mean(KKRightRateTotal);
plot(RotateAngle,KSRightRateMean,'b',RotateAngle,MSRightRateMean,'r');
end