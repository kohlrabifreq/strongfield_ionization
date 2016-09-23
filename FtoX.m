function [xpos] = FtoX(pulseLength, Aenv, p99, cosFreq)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

tNumPoints=1000;
sigEnv=p99/2/3;
mP=1.67e-27;
mSr=84*mP;
vInit=0;
xInit=0;

pulseTime=0:pulseLength/tNumPoints:pulseLength;
time=0:pulseLength/tNumPoints:pulseLength*3;
tCenter=mean(time);
waveform=@(Aenv, sigEnv, cosFreq) Aenv*exp(-(time-tCenter).^2/(2*sigEnv^2)).*cos(cosFreq.*time);

force=waveform(Aenv,sigEnv,cosFreq);
acc=force/mSr;
vel=cumtrapz(time,acc)+vInit;
xpos=cumtrapz(time,vel);

figure(1);
subplot(221);
plot(time, force);
subplot(222);
plot(time, acc);
subplot(223);
plot(time,vel);
subplot(224);
plot(time,xpos);



end

