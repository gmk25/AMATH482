%%Part 1
%% instructor given code
load handel
v = y';
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');
% plays the music
p8 = audioplayer(v,Fs);
playblocking(p8);
%% Discretizing
L = 9;
n = 73113; %chose my own n
t2 = linspace(0,L,n+1);
t = t2(1:n);
k = ((2*pi)/L)*[0:n/2 -n/2:-1]; %modified k
ks = fftshift(k); %shifting so we don't have to shift when plotting
%% Gabor filter
a = [1 8 50];
tau = 0:0.1:L;
Sgt_spec = zeros(length(tau),n);
for j = 1:length(a)
for i = 1:length(tau)
    g = exp(-a(j)*(t-tau(i)).^2); %Gabor filter
    Sg = g.*v; %applying filter to signal
    Sgt = fft(Sg); %taking FFT of filtered signal frequency
    Sgt_spec(i,:) = fftshift(abs(Sgt));
end
subplot(length(a),1,j)
colorbar
pcolor(tau,ks,Sgt_spec.')
shading interp
colormap(hot)
title(['a =',num2str(a(j))])
xlabel('time')
ylabel('frequency')
end
%% Mexican Hat Wavelet
tau = 0:3:L;
Sgt_spec = zeros(length(tau),n);
for i = 1:length(tau)
    g = (1-(t-tau(i)).^2).*exp((-(t-tau(i)).^2)/2); %Mexican Hat wavelet
    Sg = g.*v; %applying filter to signal
    Sgt = fft(Sg); %taking FFT of filtered signal frequency
    Sgt_spec(i,:) = fftshift(abs(Sgt));
end
subplot(length(a),1,j)
colorbar
pcolor(tau,ks,Sgt_spec.')
shading interp
colormap(hot)
title(['a =',num2str(a(j))])
xlabel('time')
ylabel('frequency')
%% Mother Wavelet
a = [0.1 10 100];
tau = 0:3:L;
Sgt_spec = zeros(length(tau),n);
for j = 1:length(a)
for i = 1:length(tau)
    g = 1/sqrt(a(j)).*((t-tau(i))./a(j)); %Mother wavelet function
    Sg = g.*v; %applying filter to signal
    Sgt = fft(Sg); %taking FFT of filtered signal frequency
    Sgt_spec(i,:) = fftshift(abs(Sgt));
end
subplot(length(a),1,j)
colorbar
pcolor(tau,ks,Sgt_spec.')
shading interp
colormap(hot)
title(['a =',num2str(a(j))])
xlabel('time')
ylabel('frequency')
end
%% Part 2
%% Music 1- Piano
clc;clear all;
% Instructor given code
[y,Fs] = audioread('music1.wav'); %y= sampled data (m=# of audio samples read, n=#of audio channels in file), Fs= sample rate(Hz)
tr_piano=length(y)/Fs; % record time in seconds
figure(1)
subplot(2,1,1)
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); 
ylabel('Amplitude');
title('Mary had a little lamb (piano)');
p8 = audioplayer(y,Fs);
playblocking(p8);
%% Discretizing
sig = y.';
L = tr_piano;
n = length(y); %chose my own n
t2 = linspace(0,L,n+1);
t = t2(1:n);
k = ((2*pi)/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k); %shifting so we don't have to shift when plotting
%% Creating Spectrogram
tslide = 0:0.1:L;
a1 = 100;
vector_piano = zeros();
ind_piano = zeros();
St_spec1 = zeros(length(tslide),n);
for b = 1:length(tslide)
g1 = exp(-a1*(t-tslide(b)).^2);
St1 = fft(g1.*sig);%St1 is in the frequency domain
St_spec1(b, :) = fftshift(abs(St1));
    %find max at each point
    k0 = abs(max(St1));
    indice = find(abs(St1)==k0);
    ksindexes = find(abs(ks)==k(indice(1)));
    vector_piano(:,b) = abs(ks(ksindexes(1))./(2*pi)); %convert to Hz
    ind_piano(:,b) = ksindexes(2);
end
St_spec1 = St_spec1(:,[ind_piano]);
%% Plotting
pcolor(tslide,vector_piano,St_spec1.') %unfiltered function
shading interp
colormap(hot)
colorbar
title('Piano Spectrogram')
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
set(gca,'XTick',[0:1:L])
%% Music 2- Recorder
clear;clc;
[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs; % record time in seconds
figure(2)
subplot(2,1,1)
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); 
ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
p8 = audioplayer(y,Fs);
playblocking(p8);
%% Discretizing
sig = y.';
L = tr_rec;
n = length(y); %chose my own n
t2 = linspace(0,L,n+1);
t = t2(1:n);
k = ((2*pi)/L)*[0:n/2-1 -n/2:-1];
ks = fftshift(k); %shifting so we don't have to shift when plotting
%% Creating spectrogram
tslide = 0:0.1:L;
a2 = 100;
vector_rec = zeros();
ind_rec = zeros();
St_spec2 = zeros(length(tslide),n);
for b = 1:length(tslide)
g2 = exp(-a2*(t-tslide(b)).^2);
St2 = fft(g2.*sig); %St1 is in the frequency domain
St_spec2(b, :) = fftshift(abs(St2));
    %find max at each point
    k0 = abs(max(St2));
    indice = find(abs(St2)==k0);
    ksindexes = find(abs(ks)==k(indice(1)));
    vector_rec(:,b) = abs(ks(ksindexes(1))./(2*pi)); %convert to Hz
    ind_rec(:,b) = ksindexes(2);
end
St_spec2 = St_spec2(:,[ind_rec]);
%% Plotting
pcolor(tslide,vector_rec,St_spec2.') %unfiltered function
shading interp
colormap(hot)
colorbar
title('Recorder Spectrogram')
xlabel('Time (sec)')
ylabel('Frequency (Hz)')
set(gca,'XTick',[0:1:L])
