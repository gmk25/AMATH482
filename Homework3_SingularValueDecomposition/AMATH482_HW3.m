%% Test 1
load('cam1_1.mat')
numFrames = size(vidFrames1_1,4);
dif = zeros(1,307200);
for j = 1:numFrames
X = vidFrames1_1(:,:,:,j);
X_gray = rgb2gray(X);
if j ==1
new = reshape(rgb2gray(vidFrames1_1(:,:,:,j)),[1,307200]);
end
new = imadd(reshape(rgb2gray(vidFrames1_1(:,:,:,j)),[1,307200]),new);
end
avgframe = new./numFrames;
for j = 1:numFrames
dif(j,:) = reshape(rgb2gray(vidFrames1_1(:,:,:,j)),[1,307200])-avgframe;
end
%%
X = vidFrames1_1(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows1_1 = zeros();
tempcols1_1 = zeros();
temprows1_1(1) = 253;
tempcols1_1(1) = 559;
for i = 2:numFrames
    X = vidFrames1_1(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(cols-tempcols1_1(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(rows(count)-temprows1_1(i-1,:)));
           temprows1_1(i,:) = (rows(count(rindex)));
           tempcols1_1(i,:) = (cols(count(rindex)));
end
%%
for j = 1:size(vidFrames1_1,4)
X = vidFrames1_1(:,:,:,j);
X_gray = rgb2gray(X);
imshow(X);
hold on
plot(tempcols1_1(j),temprows1_1(j),'ro')
pause(0.1)
end
%%
load('cam2_1.mat')
numFrames = size(vidFrames2_1,4);
dif = zeros(1,307200);
for j = 1:numFrames
X = vidFrames2_1(:,:,:,j);
X_gray = rgb2gray(X);
if j ==1
new = reshape(rgb2gray(vidFrames2_1(:,:,:,j)),[1,307200]);
end
new = imadd(reshape(rgb2gray(vidFrames2_1(:,:,:,j)),[1,307200]),new);
end
avgframe = new./numFrames;
for j = 1:numFrames
dif(j,:) = reshape(rgb2gray(vidFrames2_1(:,:,:,j)),[1,307200])-avgframe;
end
%%
X = vidFrames2_1(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows2_1 = zeros();
tempcols2_1 = zeros();
temprows2_1(1) = 266;
tempcols2_1(1) = 317;
for i = 2:numFrames
    X = vidFrames2_1(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(cols-tempcols2_1(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(rows(count)-temprows2_1(i-1,:)));
           temprows2_1(i,:) = (rows(count(rindex)));
           tempcols2_1(i,:) = (cols(count(rindex)));
end
%%
for j = 1:size(vidFrames2_1,4)
X = vidFrames2_1(:,:,:,j);
X_gray = rgb2gray(X);
imshow(X);
hold on
plot(tempcols2_1(j),temprows2_1(j),'ro')
pause(0.1)
end
%%
load('cam3_1.mat')
numFrames = size(vidFrames3_1,4);
dif = zeros(1,307200);
for j = 1:numFrames
X = vidFrames3_1(:,:,:,j);
X_gray = rgb2gray(X);
if j ==1
new = reshape(rgb2gray(vidFrames3_1(:,:,:,j)),[1,307200]);
end
new = imadd(reshape(rgb2gray(vidFrames3_1(:,:,:,j)),[1,307200]),new);
end
avgframe = new./numFrames;
for j = 1:numFrames
dif(j,:) = reshape(rgb2gray(vidFrames3_1(:,:,:,j)),[1,307200])-avgframe;
end
%%
X = vidFrames3_1(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows3_1 = zeros();
tempcols3_1 = zeros();
temprows3_1(1) = 271;
tempcols3_1(1) = 319;
for i = 2:numFrames
    X = vidFrames3_1(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(rows-temprows3_1(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(cols(count)-tempcols3_1(i-1,:)));
           temprows3_1(i,:) = (rows(count(rindex)));
           tempcols3_1(i,:) = (cols(count(rindex)));
end
%%
for j = 1:numFrames
X = vidFrames3_1(:,:,:,j);
X_gray = rgb2gray(X);
imshow(X);
hold on
plot(tempcols3_1(j),temprows3_1(j),'ro')
pause(0.1)
end
%% Test 2
load('cam1_2.mat')
numFrames = size(vidFrames1_2,4);
dif = zeros(1,307200);
for j = 1:numFrames
X = vidFrames1_2(:,:,:,j);
X_gray = rgb2gray(X);
if j ==1
new = reshape(rgb2gray(vidFrames1_2(:,:,:,j)),[1,307200]);
end
new = imadd(reshape(rgb2gray(vidFrames1_2(:,:,:,j)),[1,307200]),new);
end
avgframe = new./numFrames;
for j = 1:numFrames
dif(j,:) = reshape(rgb2gray(vidFrames1_2(:,:,:,j)),[1,307200])-avgframe;
end
%%
X = vidFrames1_2(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows1_2 = zeros();
tempcols1_2 = zeros();
temprows1_2(1) = 309;
tempcols1_2(1) = 326;
for i = 2:numFrames
    X = vidFrames1_2(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(cols-tempcols1_2(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(rows(count)-temprows1_2(i-1,:)));
           temprows1_2(i,:) = (rows(count(rindex)));
           tempcols1_2(i,:) = (cols(count(rindex)));
end
%%
for j = 1:size(vidFrames1_2,4)
X = vidFrames1_2(:,:,:,j);
X_gray = rgb2gray(X);
imshow(X); 
hold on
plot(tempcols1_2(j),temprows1_2(j),'ro')
pause(0.1)
end
%%
load('cam2_2.mat')
numFrames = size(vidFrames2_2,4);
dif = zeros(1,307200);
for j = 1:numFrames
X = vidFrames2_2(:,:,:,j);
X_gray = rgb2gray(X);
if j ==1
new = reshape(rgb2gray(vidFrames2_2(:,:,:,j)),[1,307200]);
end
new = imadd(reshape(rgb2gray(vidFrames2_2(:,:,:,j)),[1,307200]),new);
end
avgframe = new./numFrames;
for j = 1:numFrames
dif(j,:) = reshape(rgb2gray(vidFrames2_2(:,:,:,j)),[1,307200])-avgframe;
end
%%
X = vidFrames2_2(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows2_2 = zeros();
tempcols2_2 = zeros();
temprows2_2(1) = 361;
tempcols2_2(1) = 317;
for i = 2:numFrames
    X = vidFrames2_2(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(cols-tempcols2_2(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(rows(count)-temprows2_2(i-1,:)));
           temprows2_2(i,:) = (rows(count(rindex)));
           tempcols2_2(i,:) = (cols(count(rindex)));
end
%%
for j = 1:numFrames
X = vidFrames2_2(:,:,:,j);
X_gray = rgb2gray(X);
imshow(X); 
hold on
plot(tempcols2_2(j),temprows2_2(j),'ro')
pause(0.1)
end
%%
load('cam3_2.mat')
numFrames = size(vidFrames3_2,4);
dif = zeros(1,307200);
for j = 1:numFrames
X = vidFrames3_2(:,:,:,j);
X_gray = rgb2gray(X);
if j ==1
new = reshape(rgb2gray(vidFrames3_2(:,:,:,j)),[1,307200]);
end
new = imadd(reshape(rgb2gray(vidFrames3_2(:,:,:,j)),[1,307200]),new);
end
avgframe = new./numFrames;
for j = 1:numFrames
dif(j,:) = reshape(rgb2gray(vidFrames3_2(:,:,:,j)),[1,307200])-avgframe;
end
%%
X = vidFrames3_2(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows3_2 = zeros();
tempcols3_2 = zeros();
temprows3_2(1) = 246;
tempcols3_2(1) = 350;
for i = 2:numFrames
    X = vidFrames3_2(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(rows-temprows3_2(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(cols(count)-tempcols3_2(i-1,:)));
           temprows3_2(i,:) = (rows(count(rindex)));
           tempcols3_2(i,:) = (cols(count(rindex)));
end
%% Test 3
load('cam1_3.mat')
numFrames = size(vidFrames1_3,4);
dif = zeros(1,307200);
%find the mean
for j = 1:numFrames
X = vidFrames1_3(:,:,:,j);
X_gray = rgb2gray(X);
if j ==1
new = reshape(rgb2gray(vidFrames1_3(:,:,:,j)),[1,307200]);
end
new = imadd(reshape(rgb2gray(vidFrames1_3(:,:,:,j)),[1,307200]),new);
end
avgframe = new./numFrames;
for j = 1:numFrames
dif(j,:) = reshape(rgb2gray(vidFrames1_3(:,:,:,j)),[1,307200])-avgframe;
end
%%
X = vidFrames1_3(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows1_3 = zeros();
tempcols1_3 = zeros();
temprows1_3(1) = 306;
tempcols1_3(1) = 324;
for i = 2:numFrames
    X = vidFrames1_3(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(cols-tempcols1_3(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(rows(count)-temprows1_3(i-1,:)));
           temprows1_3(i,:) = (rows(count(rindex)));
           tempcols1_3(i,:) = (cols(count(rindex)));
end
%%
for j = 1:numFrames
X = vidFrames1_3(:,:,:,j);
X_gray = rgb2gray(X);
imshow(X);
hold on
plot(tempcols1_3(j),temprows1_3(j),'ro')
pause(0.1)
end
%%
load('cam2_3.mat')
numFrames = size(vidFrames2_3,4);
dif = zeros(1,307200);
%find the mean
for j = 1:numFrames
X = vidFrames2_3(:,:,:,j);
X_gray = rgb2gray(X);
if j ==1
new = reshape(rgb2gray(vidFrames2_3(:,:,:,j)),[1,307200]);
end
new = imadd(reshape(rgb2gray(vidFrames2_3(:,:,:,j)),[1,307200]),new);
end
avgframe = new./numFrames;
for j = 1:numFrames
dif(j,:) = reshape(rgb2gray(vidFrames2_3(:,:,:,j)),[1,307200])-avgframe;
end
%%
X = vidFrames2_3(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows2_3 = zeros();
tempcols2_3 = zeros();
temprows2_3(1) = 296;
tempcols2_3(1) = 239;
for i = 2:numFrames
    X = vidFrames2_3(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(cols-tempcols2_3(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(rows(count)-temprows2_3(i-1,:)));
           temprows2_3(i,:) = (rows(count(rindex)));
           tempcols2_3(i,:) = (cols(count(rindex)));
end
%%
for j = 1:numFrames
X = vidFrames2_3(:,:,:,j);
X_gray = rgb2gray(X);
imshow(X);
hold on
plot(tempcols2_3(j),temprows2_3(j),'ro')
pause(0.1)
end
%%
load('cam3_3.mat')
numFrames = size(vidFrames3_3,4);
dif = zeros(1,307200);
%find the mean
for j = 1:numFrames
X = vidFrames3_3(:,:,:,j);
X_gray = rgb2gray(X);
if j ==1
new = reshape(rgb2gray(vidFrames3_3(:,:,:,j)),[1,307200]);
end
new = imadd(reshape(rgb2gray(vidFrames3_3(:,:,:,j)),[1,307200]),new);
end
avgframe = new./numFrames;
for j = 1:numFrames
dif(j,:) = reshape(rgb2gray(vidFrames3_3(:,:,:,j)),[1,307200])-avgframe;
end
%%
X = vidFrames3_3(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows3_3 = zeros();
tempcols3_3 = zeros();
temprows3_3(1) = 233;
tempcols3_3(1) = 357;
for i = 2:numFrames
    X = vidFrames3_3(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(rows-temprows3_3(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(cols(count)-tempcols3_3(i-1,:)));
           temprows3_3(i,:) = (rows(count(rindex)));
           tempcols3_3(i,:) = (cols(count(rindex)));
end
%% Test 4
load('cam1_4.mat')
numFrames = size(vidFrames1_4,4);
dif = zeros(1,307200);
for j = 1:numFrames
X = vidFrames1_4(:,:,:,j);
X_gray = rgb2gray(X);
if j ~=1
dif(j,:) = reshape(rgb2gray(vidFrames1_4(:,:,:,j)) - rgb2gray(vidFrames1_4(:,:,:,j-1)),[1,307200]);
end
end
%%
X = vidFrames1_4(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows1_4 = zeros();
tempcols1_4 = zeros();
temprows1_4(1) = 305;
tempcols1_4(1) = 376;
for i = 2:numFrames
    X = vidFrames1_4(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(cols-tempcols1_4(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(rows(count)-temprows1_4(i-1,:)));
           [d, cindex] = min(abs(cols(count)-tempcols1_4(i-1,:)));
           temprows1_4(i,:) = (rows(count(rindex)));
           tempcols1_4(i,:) = (cols(count(cindex)));
end
%%
for j = 1:numFrames
X = vidFrames1_4(:,:,:,j);
X_gray = rgb2gray(X);
imshow(X);
hold on
plot(tempcols1_4(j),temprows1_4(j),'ro')
pause(0.1)
end
%%
load('cam2_4.mat')
numFrames = size(vidFrames2_4,4);
dif = zeros(1,307200);
for j = 1:numFrames
X = vidFrames2_4(:,:,:,j);
X_gray = rgb2gray(X);
if j ~=1
dif(j,:) = reshape(rgb2gray(vidFrames2_4(:,:,:,j)) - rgb2gray(vidFrames2_4(:,:,:,j-1)),[1,307200]);
end
end
%%
X = vidFrames2_4(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows2_4 = zeros();
tempcols2_4 = zeros();
temprows2_4(1) = 244;
tempcols2_4(1) = 244;
for i = 2:numFrames
    X = vidFrames2_4(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(cols-tempcols2_4(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(rows(count)-temprows2_4(i-1,:)));
           [d, cindex] = min(abs(cols(count)-tempcols2_4(i-1,:)));
           temprows2_4(i,:) = (rows(count(rindex)));
           tempcols2_4(i,:) = (cols(count(cindex)));
end
%%
for j = 1:numFrames
X = vidFrames2_4(:,:,:,j);
X_gray = rgb2gray(X);
imshow(X);
hold on
plot(tempcols2_4(j),temprows2_4(j),'ro')
pause(0.1)
end
%%
load('cam3_4.mat')
numFrames = size(vidFrames3_4,4);
dif = zeros(1,307200);
for j = 1:numFrames
X = vidFrames3_4(:,:,:,j);
X_gray = rgb2gray(X);
if j ~=1
dif(j,:) = reshape(rgb2gray(vidFrames3_4(:,:,:,j)) - rgb2gray(vidFrames3_4(:,:,:,j-1)),[1,307200]);
end
end
%%
X = vidFrames3_4(:,:,:,1);
X_gray = rgb2gray(X);
imshow(X_gray)
ginput(1)
%% find greatest difference between frames
temprows3_4 = zeros();
tempcols3_4 = zeros();
temprows3_4(1) = 208;
tempcols3_4(1) = 399;
for i = 2:numFrames
    X = vidFrames3_4(:,:,:,i);
    X_gray = rgb2gray(X);
    temp = find(dif(i,:) == max(abs(dif(i,:))));
        [rows,cols] = find(X_gray==X_gray(temp(1)));
           count1 = abs(rows-temprows3_4(i-1));
           count = find(count1 == min(count1));
           [c, rindex] = min(abs(rows(count)-temprows3_4(i-1,:)));
           [d, cindex] = min(abs(cols(count)-tempcols3_4(i-1,:)));
           temprows3_4(i,:) = (rows(count(rindex)));
           tempcols3_4(i,:) = (cols(count(cindex)));
end
%%
for j = 1:numFrames
X = vidFrames3_4(:,:,:,j);
X_gray = rgb2gray(X);
imshow(X);
hold on
plot(tempcols3_4(j),temprows3_4(j),'ro')
pause(0.1)
end
%% Creating snapshot matrix and peforming the SVD on all tests
%% Test 1
figure(1) 
plot(1:length(temprows1_1),temprows1_1,'k-')
figure(2) 
plot(1:length(temprows2_1),temprows2_1,'k-')
figure(3) 
plot( 1:length(tempcols3_1),tempcols3_1,'k-')
%120 points: camera 1, x=11:131 camera 2: x = 19:139 camera 3: x = 9:129
parsed1_1 = [tempcols1_1(13:133) temprows1_1(13:133)];
parsed2_1 = [tempcols2_1(19:139) temprows2_1(19:139)];
parsed3_1 = [tempcols3_1(9:129) temprows3_1(9:129)];
snapshot_1 = [parsed1_1.';parsed2_1.';parsed3_1.'];
%subtracting the mean from each row
snapshot_1 = snapshot_1 - [mean(snapshot_1(1,:));mean(snapshot_1(2,:));mean(snapshot_1(3,:));mean(snapshot_1(4,:));mean(snapshot_1(5,:));mean(snapshot_1(6,:))];
[U_1,S_1,V_1] = svd(snapshot_1,'econ');
sig1 = diag(S_1);
%% Test 2
figure(1) 
plot(1:length(temprows1_2),temprows1_2,'k-')
figure(2) 
plot(1:length(temprows2_2),temprows2_2,'k-')
figure(3) 
plot(1:length(tempcols3_2),tempcols3_2,'k-')
% 71 points: camera 1, x=1:71 camera 2: x = 1:71 camera 3: x = 1:71
parsed1_2 = [tempcols1_2(1:71) temprows1_2(1:71)];
parsed2_2 = [tempcols2_2(1:71) temprows2_2(1:71)];
parsed3_2 = [tempcols3_2(1:71) temprows3_2(1:71)];
snapshot_2 = [parsed1_2.';parsed2_2.';parsed3_2.'];
snapshot_2= snapshot_2 - [mean(snapshot_2(1,:));mean(snapshot_2(2,:));mean(snapshot_2(3,:));mean(snapshot_2(4,:));mean(snapshot_2(5,:));mean(snapshot_2(6,:))];
[U_2,S_2,V_2] = svd(snapshot_2,'econ');
sig2 = diag(S_2);
%% Test 3
figure(1) 
plot(1:length(temprows1_3),temprows1_3,'k-')
figure(2)
plot(1:length(temprows2_3),temprows2_3,'k-')
figure(3)
plot( 1:length(tempcols3_3),tempcols3_3,'k-')
%104 points: 1:104
parsed1_3 = [tempcols1_3(1:21) temprows1_3(1:21)];
parsed2_3 = [tempcols2_3(1:21) temprows2_3(1:21)];
parsed3_3 = [tempcols3_3(1:21) temprows3_3(1:21)];
snapshot_3 = [parsed1_3.';parsed2_3.';parsed3_3.'];
snapshot_3 = snapshot_3 - [mean(snapshot_3(1,:));mean(snapshot_3(2,:));mean(snapshot_3(3,:));mean(snapshot_3(4,:));mean(snapshot_3(5,:));mean(snapshot_3(6,:))];
[U_3,S_3,V_3] = svd(snapshot_3,'econ');
sig3 = diag(S_3);
%% Test 4
figure(1)
plot(1:length(temprows1_4),temprows1_4,'k-')
figure(2)
plot(1:length(temprows2_4),temprows2_4,'k-')
figure(3)
plot( 1:length(tempcols3_4),tempcols3_4,'k-')
%48 points: 23:71
parsed1_4 = [tempcols1_4(6:54) temprows1_4(6:54)];
parsed2_4 = [tempcols2_4(7:55) temprows2_4(7:55)];
parsed3_4 = [tempcols3_4(23:71) temprows3_4(23:71)];
snapshot_4 = [parsed1_4.';parsed2_4.';parsed3_4.'];
snapshot_4 = snapshot_4 - [mean(snapshot_4(1,:));mean(snapshot_4(2,:));mean(snapshot_4(3,:));mean(snapshot_4(4,:));mean(snapshot_4(5,:));mean(snapshot_4(6,:))];
[U_4,S_4,V_4] = svd(snapshot_4,'econ');
sig4 = diag(S_4);
%% Plotting energies
subplot(2,2,1)
semilogy(sig1.^2/sum(sig1.^2),'ro')
title('Test 1 Singular Value Energies')
xlabel('Mode')
ylabel('Energy')
subplot(2,2,2)
semilogy(sig2.^2/sum(sig2.^2),'ro')
title('Test 2 Singular Value Energies')
xlabel('Mode')
ylabel('Energy')
subplot(2,2,3)
semilogy(sig3.^2/sum(sig3.^2),'ro')
title('Test 3 Singular Value Energies')
xlabel('Mode')
ylabel('Energy')
subplot(2,2,4)
semilogy(sig4.^2/sum(sig4.^2),'ro')
title('Test 4 Singular Value Energies')
xlabel('Mode')
ylabel('Energy')
%% plotting x- and y-coordinates
subplot(2,2,1)
plot(1:length(temprows1_1),temprows1_1,'r')
hold on
plot(1:length(temprows2_1),temprows2_1,'b')
plot(1:length(temprows3_1),temprows3_1,'g')
title('Test 1')
xlabel('Time')
ylabel('Vertical Position')
legend('Camera 1', 'Camera 2', 'Camera 3')
subplot(2,2,2)
plot(1:length(temprows1_2),temprows1_2,'r')
hold on
plot(1:length(temprows2_2),temprows2_2,'b')
plot(1:length(temprows3_2),temprows3_2,'g')
title('Test 2')
xlabel('Time')
ylabel('Vertical Position')
legend('Camera 1', 'Camera 2', 'Camera 3')
subplot(2,2,3)
plot(1:length(temprows1_3),temprows1_3,'r')
hold on
plot(1:length(temprows2_3),temprows2_3,'b')
plot(1:length(temprows3_3),temprows3_3,'g')
title('Test 3')
xlabel('Time')
ylabel('Vertical Position')
legend('Camera 1', 'Camera 2', 'Camera 3')
subplot(2,2,4)
plot(1:length(temprows1_4),temprows1_4,'r')
hold on
plot(1:length(temprows2_4),temprows2_4,'b')
plot(1:length(temprows3_4),temprows3_4,'g')
title('Test 4')
xlabel('Time')
ylabel('Vertical Position')
legend('Camera 1', 'Camera 2', 'Camera 3')
%%
subplot(2,2,1)
plot(1:length(tempcols1_1),tempcols1_1,'r')
hold on
plot(1:length(tempcols2_1),tempcols2_1,'b')
plot(1:length(tempcols3_1),tempcols3_1,'g')
title('Test 1')
xlabel('Time')
ylabel('Horizontal Position')
legend('Camera 1', 'Camera 2', 'Camera 3')
subplot(2,2,2)
plot(1:length(tempcols1_2),tempcols1_2,'r')
hold on
plot(1:length(tempcols2_2),tempcols2_2,'b')
plot(1:length(tempcols3_2),tempcols3_2,'g')
title('Test 2')
xlabel('Time')
ylabel('Horizontal Position')
legend('Camera 1', 'Camera 2', 'Camera 3')
subplot(2,2,3)
plot(1:length(tempcols1_3),tempcols1_3,'r')
hold on
plot(1:length(tempcols2_3),tempcols2_3,'b')
plot(1:length(tempcols3_3),tempcols3_3,'g')
title('Test 3')
xlabel('Time')
ylabel('Horizontal Position')
legend('Camera 1', 'Camera 2', 'Camera 3')
subplot(2,2,4)
plot(1:length(tempcols1_4),tempcols1_4,'r')
hold on
plot(1:length(tempcols2_4),tempcols2_4,'b')
plot(1:length(tempcols3_4),tempcols3_4,'g')
title('Test 4')
xlabel('Time')
ylabel('Horizontal Position')
legend('Camera 1', 'Camera 2', 'Camera 3')
