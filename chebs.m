%AUTHOR:Ashish Rai
%(DESIGN AND IMPLEMENT LOW PASS BUTTERWORTH
%FILTER USING THE GIVEN PARAMETERS)
%DATE:12/10/2015
home;
rp=0.3;%pass band ripple
rs=60;%stop band ripple
wp=1000;%pass band frequency
ws=2000;%stop band frequency
fs=9000;%sampling frequency
w1=2*wp/fs;
w2=2*ws/fs; % for digital to analog conversion
[n,wn]=buttord(w1,w2,rp,rs,'s'); %'s' for analog filter
[z,p,k]=butter(n,wn);
[b,a] = butter(n,wn,'s');
[h,w]=freqz(a,b,512,1000); %converted to analog
plot(w/pi,20*log10(abs(h)),'linewidth',2);
xlabel('Frequency-->');
ylabel('Analog Frequency response');
title('Butterworth Lowpass Filter')
grid on;
[numz,denz]=bilinear(b,a,0.5);%transfrom to digital



