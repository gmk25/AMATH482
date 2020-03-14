%loading test data file
clear; close all; clc;
load('Testdata.mat')
%% Starter code given by instructor
L=15; % defining spatial domain
n=64; % setting Fourier modes
x2=linspace(-L,L,n+1);
x=x2(1:n); %consider only the first n points
y=x; %y-axis will be defined the same way as the x-axis
z=x;%z-axis will be defined the same way as the x-axis
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; %k rescaled to 2pi domain
ks=fftshift(k); %shift Fourier transform
[X,Y,Z]=meshgrid(x,y,z); %x, y, and z axis
[Kx,Ky,Kz]=meshgrid(ks,ks,ks); %fourier transform shift of x, y and z axes on Fourier transform
for j=1:20
Un(:,:,:)=reshape(Undata(j,:),n,n,n); %reshape Undata(j,:) into a nxnxn matrix
close all, isosurface(X,Y,Z,abs(Un),0.4)
axis([-20 20 -20 20 -20 20]), grid on, drawnow
pause(1)
end
xlabel('x position')
ylabel('y position')
zlabel('z position')
set(gca,'Fontsize',16)
%%
test = zeros(64,64,64);
for jj = 1:20
    Un(:,:,:)=reshape(Undata(jj,:),n,n,n); %reshape Undata(j,:) into a nxnxn matrix
    test = test+Un;
end
%% Averaging
ut = fftn(test); %taking FFT of the sum of all Un, that is 20 points
ave = abs(ut)/20; %averaging the fftshift of ut
isosurface(Kx,Ky,Kz,fftshift(ave)/fftshift(max(ave,[],'all')),.99)%Plotting frequencies while scaling the average height.
xlabel('x frequency')
ylabel('y frequency')
zlabel('z frequency')
set(gca,'Fontsize',16)
grid on, drawnow
%Finding the central frequency
maxave = max(ave,[],'all');
for x1 = 1:64
    for y1 = 1:64
        for z1 = 1:64
            if ave(x1,y1,z1)==maxave
                k1 = x1;
                k2 = y1;
                k3 = z1;
            end
        end
    end
end
%% Filtering
tau = 2;
gfil  = exp(-tau*(((Kx-Kx(k1,k2,k3)).^2))) + exp(-tau*(((Ky-Ky(k1,k2,k3)).^2))) + exp(-tau*(((Kz-Kz(k1,k2,k3)).^2)));
%Plotting the marble's positiion at each observation
Xvec = zeros(20,1);
Yvec = zeros(20,1);
Zvec = zeros(20,1);
for i = 1:20
Un(:,:,:)=reshape(Undata(i,:),n,n,n);
uft = gfil.*(fftn(Un)); %Apply filter to the signal in frequency space ifftshift(ave)
uf = (ifftn(uft)); %Convert filtered signal back to time domain
isosurface(X,Y,Z,abs(uf),1)
[f,v] = isosurface(X,Y,Z,abs(uf),1); %Plotting path of marble against filtered signal
Xvec(i,1) = v(1,1);
Yvec(i,1) = v(1,2);
Zvec(i,1) = v(1,3);
grid on, drawnow
hold on
end
xlabel('x position')
ylabel('y position')
zlabel('z position')
set(gca,'Fontsize',16)
%% Plotting the path
plot3(Xvec,Yvec,Zvec)
grid on, drawnow
xlabel('x position')
ylabel('y position')
zlabel('z position')
set(gca,'Fontsize',16)
