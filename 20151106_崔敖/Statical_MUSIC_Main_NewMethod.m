clear;
%clc;
%% Init and Settings
RotateAngle=0;

ExpSetting=struct('SampleFrequency',48000, ...
                    'WaveSpeed',340, ...
                    'SNR',10, ...
                    'TimeSnaps',5000);

SignalSetting=struct( 'SignalNum',2, ...
                        'SignalFrequency',[1000,1000], ...
                        'IncidentAngle',[round(unifrnd(0,180)),round(unifrnd(0,180))], ...
                        'InitialPhase',[0,pi], ...
                        'Power',[0,0]);

% SignalSetting=struct( 'SignalNum',2, ...
%                         'SignalFrequency',[1000,1000], ...
%                         'IncidentAngle',[45,55], ...
%                         'InitialPhase',[0,pi], ...
%                         'Power',[0,0]);
                    
% s=round(unifrnd(0,180));
% SignalSetting=struct( 'SignalNum',2, ...
%                         'SignalFrequency',[1000,1000], ...
%                         'IncidentAngle',[s,s], ...
%                         'InitialPhase',[0,0], ...
%                         'Power',[0,0]);
                    
                    
ArraySetting=struct(  'ElementNum',3, ...
                        'ElementSpacingTime',2, ...
                        'AssistArrayElementSpacingTime',2,....
                        'AssistArrayRotateAngle',0);

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
                        'AssistSteeringVectors_N',double.empty(), ...
                        'SpectrumAngles',0:0.5:180, ...
                        'MusicSpectrum',double.empty(), ...
                        'AssistSpectrumAngles',0:0.5:180, ...
                        'AssistMusicSpectrum',double.empty(), ...
                        'AssistFixedMusicSpectrum',double.empty());

Array=struct(         'ElementNum',ArraySetting.ElementNum, ...
                        'ElementSpacing',ArraySetting.ElementSpacingTime*min(Signal.WaveLength), ...
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
                        'MusicBound',double.empty());

AssistArray_P=struct(   'ElementNum',ArraySetting.ElementNum, ...
                        'ElementSpacing',ArraySetting.AssistArrayElementSpacingTime*min(Signal.WaveLength), ...
                        'SteeringAngles',0:0.5:180, ...
                        'SteeringSpacing',double.empty(), ...
                        'SteeringVectors',double.empty(), ...
                        'MusicSpectrum',double.empty(), ...
                        'FixedMusicSpectrum',double.empty(), ...
                        'FixedMusicSpectrumMean',double.empty(), ...
                        'FixedMusicVar',double.empty(), ...
                        'FixedMusicBound',double.empty(), ...
                        'RxSignal', double.empty(), ...
                        'RMatrixTime',double.empty(), ...
                        'RMatrix',double.empty(), ...
                        'NoiseBase',double.empty(), ...
                        'MusicSpectrumMean',double.empty(), ...
                        'MusicVar',double.empty(), ...
                        'MusicBound',double.empty());

AssistArray_N=struct(   'ElementNum',ArraySetting.ElementNum, ...
                        'ElementSpacing',ArraySetting.AssistArrayElementSpacingTime*min(Signal.WaveLength), ...
                        'SteeringAngles',0:0.5:180, ...
                        'SteeringSpacing',double.empty(), ...
                        'SteeringVectors',double.empty(), ...
                        'MusicSpectrum',double.empty(), ...
                        'FixedMusicSpectrum',double.empty(), ...
                        'FixedMusicSpectrumMean',double.empty(), ...
                        'FixedMusicVar',double.empty(), ...
                        'FixedMusicBound',double.empty(), ...
                        'RxSignal', double.empty(), ...
                        'RMatrixTime',double.empty(), ...
                        'RMatrix',double.empty(), ...
                        'NoiseBase',double.empty(), ...
                        'MusicSpectrumMean',double.empty(), ...
                        'MusicVar',double.empty(), ...
                        'MusicBound',double.empty());
tic
while RotateAngle<=180
  ArraySetting.AssistArrayRotateAngle=RotateAngle;
  Signal.RotateIncidentAngle_P=SignalSetting.IncidentAngle+ArraySetting.AssistArrayRotateAngle;
  Signal.RotateIncidentAngle_N=SignalSetting.IncidentAngle-ArraySetting.AssistArrayRotateAngle;
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

  Signal.SteeringSpacing=2*pi*SignalSetting.SignalFrequency(1).*Array.ElementSpacing.*cosd(Signal.IncidentAngle)./ExpSetting.WaveSpeed;    for Counter=1:1:Array.ElementNum
      Signal.SteeringVectors(Counter,:)=exp(-1i*(Counter-1).*Signal.SteeringSpacing);
  end

  Signal.AssistSteeringSpacing_P=2*pi*SignalSetting.SignalFrequency(1).*AssistArray_P.ElementSpacing.*cosd(Signal.RotateIncidentAngle_P)./ExpSetting.WaveSpeed;
  for Counter=1:1:Array.ElementNum
    Signal.AssistSteeringVectors_P(Counter,:)=exp(-1i*(Counter-1).*Signal.AssistSteeringSpacing_P);
  end

  Signal.AssistSteeringSpacing_N=2*pi*SignalSetting.SignalFrequency(1).*AssistArray_N.ElementSpacing.*cosd(Signal.RotateIncidentAngle_N)./ExpSetting.WaveSpeed;
  for Counter=1:1:Array.ElementNum
    Signal.AssistSteeringVectors_N(Counter,:)=exp(-1i*(Counter-1).*Signal.AssistSteeringSpacing_N);
  end

  for Counter=1:1:ExpSetting.TimeSnaps
    Array.RxSignal(:,Counter)=Signal.SteeringVectors*sqrt(2)*[cos(Signal.DigitalFrequency(1)*Counter);sin(Signal.DigitalFrequency(1)*Counter)];
%     Array.RxSignal(:,Counter)=Signal.SteeringVectors*sqrt(2)*[cos(Signal.DigitalFrequency(1)*Counter)];

  end

  for Counter=1:1:Array.ElementNum
    Array.RxSignal(Counter,:)=awgn(Array.RxSignal(Counter,:),ExpSetting.SNR,'measured');
  end

  for Counter=1:1:ExpSetting.TimeSnaps
    AssistArray_P.RxSignal(:,Counter)=Signal.AssistSteeringVectors_P*sqrt(2)*[cos(Signal.DigitalFrequency(1)*Counter);sin(Signal.DigitalFrequency(1)*Counter)];
%     AssistArray_P.RxSignal(:,Counter)=Signal.AssistSteeringVectors_P*sqrt(2)*[cos(Signal.DigitalFrequency(1)*Counter)];

  end
  for Counter=1:1:AssistArray_P.ElementNum
    AssistArray_P.RxSignal(Counter,:)=awgn(AssistArray_P.RxSignal(Counter,:),ExpSetting.SNR,'measured');
  end

  for Counter=1:1:ExpSetting.TimeSnaps
    AssistArray_N.RxSignal(:,Counter)=Signal.AssistSteeringVectors_N*sqrt(2)*[cos(Signal.DigitalFrequency(1)*Counter);sin(Signal.DigitalFrequency(1)*Counter)];
%     AssistArray_N.RxSignal(:,Counter)=Signal.AssistSteeringVectors_N*sqrt(2)*[cos(Signal.DigitalFrequency(1)*Counter)];

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
    Array.MusicSpectrum(RotateAngle+1,Counter)=abs(inv((Array.SteeringVectors(:,Counter))'*Array.NoiseBase*Array.NoiseBase'*Array.SteeringVectors(:,Counter)));
  end
  Array.MusicSpectrum(RotateAngle+1,:)=Array.MusicSpectrum(RotateAngle+1,:)./max(Array.MusicSpectrum(RotateAngle+1,:));
  
  for Counter=1:1:size(AssistArray_P.SteeringAngles,2)
    AssistArray_P.MusicSpectrum(RotateAngle+1,Counter)=abs(inv((AssistArray_P.SteeringVectors(:,Counter))'*AssistArray_P.NoiseBase*AssistArray_P.NoiseBase'*AssistArray_P.SteeringVectors(:,Counter)));
  end

  for Counter=1:1:size(AssistArray_N.SteeringAngles,2)
    AssistArray_N.MusicSpectrum(RotateAngle+1,Counter)=abs(inv((AssistArray_N.SteeringVectors(:,Counter))'*AssistArray_N.NoiseBase*AssistArray_N.NoiseBase'*AssistArray_N.SteeringVectors(:,Counter)));
  end
  AssistArray_P.MusicSpectrum(RotateAngle+1,:)=AssistArray_P.MusicSpectrum(RotateAngle+1,:)./max(AssistArray_P.MusicSpectrum(RotateAngle+1,:));
  AssistArray_P.FixedMusicSpectrum(RotateAngle+1,:)=Array.MusicSpectrum(RotateAngle+1,:).*AssistArray_P.MusicSpectrum(RotateAngle+1,:);

  AssistArray_N.MusicSpectrum(RotateAngle+1,:)=AssistArray_N.MusicSpectrum(RotateAngle+1,:)./max(AssistArray_N.MusicSpectrum(RotateAngle+1,:));
  
  AssistArray_N.FixedMusicSpectrum(RotateAngle+1,:)=(AssistArray_P.MusicSpectrum(RotateAngle+1,:).*AssistArray_N.MusicSpectrum(RotateAngle+1,:)).^(1/2);
  AssistArray_N.FixedMusicSpectrum(RotateAngle+1,:)=AssistArray_N.FixedMusicSpectrum(RotateAngle+1,:)./max(AssistArray_N.FixedMusicSpectrum(RotateAngle+1,:));
  
  
  
%   ResultSpectrum(RotateAngle+1,:)=(AssistArray_P.MusicSpectrum(RotateAngle+1,:).*AssistArray_N.MusicSpectrum(RotateAngle+1,:).*Array.MusicSpectrum(RotateAngle+1,:)).^(1/3);
  ResultSpectrum(RotateAngle+1,:)=(AssistArray_P.MusicSpectrum(RotateAngle+1,:).*AssistArray_N.MusicSpectrum(RotateAngle+1,:)).^(1/2);
  ResultSpectrum(RotateAngle+1,:)=ResultSpectrum(RotateAngle+1,:)./max(ResultSpectrum(RotateAngle+1,:));
  RotateAngle=RotateAngle+1;
end
toc
Res_P=10*log10(AssistArray_P.FixedMusicSpectrum);
Res_N=10*log10(AssistArray_N.FixedMusicSpectrum);
Res_F=10*log10(ResultSpectrum);
Res_O=10*log10(Array.MusicSpectrum);
% ST_Sp=[ResultSpectrum(:,1:Signal.IncidentAngle(1)/0.5),ResultSpectrum(:,Signal.IncidentAngle(1)/0.5+2:Signal.IncidentAngle(2)/0.5),ResultSpectrum(:,Signal.IncidentAngle(2)/0.5+2:end)];
% ST_Sp=[ResultSpectrum(:,1:Signal.IncidentAngle(1)/0.5-10),ResultSpectrum(:,Signal.IncidentAngle(1)/0.5+12:Signal.IncidentAngle(2)/0.5-10),ResultSpectrum(:,Signal.IncidentAngle(2)/0.5+12:end)];
% 
% [~,MinIndex]=min(sum(ST_Sp,2));

% figure(1)
% contour(0:0.5:180,0:180,Res_P);
% set(gca,'XTick',0:10:180,'YTick',0:10:180);
% xlabel('Rotate Angle /Deg.'), ylabel('Incident Angle /Deg.');

% figure(1)
% contour(0:0.5:180,0:180,10*log10(AssistArray_P.MusicSpectrum));
% set(gca,'XTick',0:10:180,'YTick',0:10:180);
% xlabel('Incident Angle /Deg.'), ylabel('Rotate Angle /Deg.');
% 
% figure(2)
% contour(0:0.5:180,0:180,10*log10(AssistArray_N.MusicSpectrum));
% set(gca,'XTick',0:10:180,'YTick',0:10:180);
% xlabel('Incident Angle /Deg.'), ylabel('Rotate Angle /Deg.');
% 
% figure(3)
% contour(0:0.5:180,0:180,10*log10(AssistArray_N.FixedMusicSpectrum));
% set(gca,'XTick',0:10:180,'XLim',[0,180],'YTick',0:10:180,'YLim',[0,180]);
% xlabel('Incident Angle /Deg.'), ylabel('Rotate Angle /Deg.');
% 
% figure(4)
% contour(0:0.5:180,0:180,10*log10(ResultSpectrum));
% set(gca,'XTick',0:10:180,'XLim',[0,180],'YTick',0:10:180,'YLim',[0,180]);
% xlabel('Incident Angle /Deg.'), ylabel('Rotate Angle /Deg.');

% figure(2)
% contour(0:0.5:180,0:180,Res_F);
% set(gca,'XTick',0:10:180,'YTick',0:10:180);
% xlabel('Incident Angle /Deg.'), ylabel('Rotate Angle /Deg.');

% figure(3)
% plot(0:0.5:180,Res_F(46,:),0:0.5:180,10*log10(Array.MusicSpectrum(46,:)));
% set(gca,'XTick',0:10:180,'XLim',[0,180],'YTick',-120:10:0);
% xlabel('Incident Angle /Deg.'), ylabel('MUSIC Spectrum /dB');

% figure(5)
% plot(0:0.5:180,AssistArray_N.MusicSpectrum(46,:),'r:',0:0.5:180,AssistArray_P.MusicSpectrum(46,:),'g:',0:0.5:180,Array.MusicSpectrum(46,:),'b:');
% set(gca,'XTick',0:5:180,'XLim',[0,180],'YTick',0:0.02:1,'YLim',[0,1]);
% xlabel('Incident Angle /Deg.'), ylabel('Improved and Original MUSIC Spectrum @ Rotate Angle=46 Deg. /dB');

Index=ceil(acosd(0.5/2/ArraySetting.ElementSpacingTime))+1;
% figure(5)
% plot(0:0.5:180,10*log10(ResultSpectrum(Index,:)),'r',0:0.5:180,10*log10(AssistArray_N.FixedMusicSpectrum(Index,:)),'r:',0:0.5:180,10*log10(ResultSpectrum(91,:)),'g:',0:0.5:180,10*log10(Array.MusicSpectrum(Index,:)),'b:');
% set(gca,'XTick',0:5:180,'XLim',[0,180],'YTick',-60:10:0,'YLim',[-60,0]);
% xlabel('Incident Angle /Deg.'), ylabel('Improved and Original MUSIC Spectrum @ Rotate Angle=46 Deg. /dB');

figure(1)
plot(0:0.5:180,10*log10(ResultSpectrum(Index+1,:)),'r',0:0.5:180,10*log10(Array.MusicSpectrum(Index+1,:)),'b:');
set(gca,'XTick',0:5:180,'XLim',[0,180],'YTick',-60:10:0,'YLim',[-60,0]);
pleg1=legend('Improved DOA','Original DOA');
title('Improved and Original MUSIC Spectrum @ Rotate Angle=84');
xlabel('Incident Angle /Deg.'), ylabel('MUSIC Spetrum /dB');

figure(2)
plot(0:0.5:180,10*log10(ResultSpectrum(21,:)),'r',0:0.5:180,10*log10(Array.MusicSpectrum(21,:)),'b:');
set(gca,'XTick',0:5:180,'XLim',[0,180],'YTick',-60:10:0,'YLim',[-60,0]);
pleg2=legend('Improved DOA','Original DOA');
title('Improved and Original MUSIC Spectrum @ Rotate Angle=20');
xlabel('Incident Angle /Deg.'), ylabel('MUSIC Spetrum /dB');

figure(3)
plot(0:0.5:180,10*log10(ResultSpectrum(41,:)),'r',0:0.5:180,10*log10(Array.MusicSpectrum(41,:)),'b:');
set(gca,'XTick',0:5:180,'XLim',[0,180],'YTick',-60:10:0,'YLim',[-60,0]);
pleg3=legend('Improved DOA','Original DOA');
title('Improved and Original MUSIC Spectrum @ Rotate Angle=40');
xlabel('Incident Angle /Deg.'), ylabel('MUSIC Spetrum /dB');

figure(4)
plot(0:0.5:180,10*log10(ResultSpectrum(71,:)),'r',0:0.5:180,10*log10(Array.MusicSpectrum(71,:)),'b:');
set(gca,'XTick',0:5:180,'XLim',[0,180],'YTick',-60:10:0,'YLim',[-60,0]);
pleg4=legend('Improved DOA','Original DOA');
title('Improved and Original MUSIC Spectrum @ Rotate Angle=70');
xlabel('Incident Angle /Deg.'), ylabel('MUSIC Spetrum /dB');
% figure(6)
% plot(0:1:180,sum(ResultSpectrum(1:181,:),2));
% set(gca,'XTick',0:5:180,'XLim',[0,180],'YTick',1:0.5:5,'YLim',[1,5]);
% xlabel('Rotate Angle /Deg.'), ylabel('Sum of MUSIC Spectrum');

disp('______________________________________________');
disp('Real Incident Angles:');
disp(SignalSetting.IncidentAngle(1));
disp(SignalSetting.IncidentAngle(2));

disp('______________________________________________');
disp('Incident Angles in Orignal MUSIC Spectrum:');
CheckOrg=Array.MusicSpectrum(Index,:);
[~,Peak1]=max(CheckOrg);
CheckOrg(Peak1)=CheckOrg(Peak1)-CheckOrg(Peak1);
[~,Peak2]=max(CheckOrg);
disp((Peak1-1)*0.5);
disp((Peak2-1)*0.5);

disp('______________________________________________');
disp('Incident Angles in Improved MUSIC Spectrum:');
CheckImp=ResultSpectrum(Index,:);
[~,Peak1]=max(CheckImp);
CheckImp(Peak1)=CheckImp(Peak1)-CheckImp(Peak1);
[~,Peak2]=max(CheckImp);
disp((Peak1-1)*0.5);
disp((Peak2-1)*0.5);
% figure(4)
% plot(0:0.5:180,10*log10(Array.MusicSpectrum(46,:)));
% set(gca,'XTick',0:10:180,'YTick',-120:10:0);
% xlabel('Incident Angle /Deg.'), ylabel('MUSIC Spectrum /dB');