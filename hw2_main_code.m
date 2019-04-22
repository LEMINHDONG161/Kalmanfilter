%
% Point target detection edited by Jinho Bae
% 2008. 10. 5
% 
% ----------------------------------------------------

% Open 2-D data file
%% Home work 2
close all
clear all
load target_hw2.mat

imagesc(mf_im);
colormap(gray)
S=[0.1 0.15 0.1;0.15 0.2 0.05; 0.1 0.15 0.1]
V1=[1 0.2 0.04;0.2 1 0.2; 0.04 0.2 1]
V2=[0.1 0.02 0.004;0.02 0.1 0.02; 0.004 0.02 0.1]
V3=[0.01 0.002 0.0004;0.002 0.01 0.002; 0.0004 0.002 0.01]
K1=inv(V1);
H1=K1*S;
Y1=filter2(H1,mf_im);
K2=inv(V2);
H2=K2*S;
Y2=filter2(H2,mf_im);
K3=inv(V3);
H3=K3*S;
Y3=filter2(H3,mf_im);

M1=imbinarize(Y1,0.8);
figure,
imagesc(M1);
colormap(gray)
M1=imbinarize(Y1,0.99);
figure,
imagesc(M1);
colormap(gray)

M1=imbinarize(Y1,1.3);
figure,
imagesc(M1);
colormap(gray)

M1=imbinarize(Y1,1.4);
figure,
imagesc(M1);
colormap(gray)

M1=imbinarize(Y1,1.5);
figure,
imagesc(M1);
colormap(gray)




