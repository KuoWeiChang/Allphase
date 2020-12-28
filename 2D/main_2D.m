%test 2D mix 
clear all;
close all;
clc;
pkg load signal;
N=512;
nx=[0:N-1];
ny=[0:N-1]';
fx=[31.58,69.12,214.77];
fy=[89.61,11.52,43.98];
A=[1,1,1];
phi=[0,0.5,1.2];
x=exp_maker2(fx,fy,A,phi,N,nx,ny);
x=x+(randn(size(x))*0.5+randn(size(x))*0.5*1i)/sqrt(2);
h=hann(N);
method='polyAMI2';
[efx,efy,eA,ephi]=AMI_2D_all(x,N,nx,ny,method,h);
effinal=zeros(length(efx),2);
for k=1:length(efx)
  efxk=efx(k);
  efyk=efy(k);
  for iter=1:5
    [efxk,efyk]=polyAMI2D(x,N,nx,ny,efxk,efyk,1);
    [efxk,efyk]=polyAMI2D(x,N,nx,ny,efxk,efyk,2);
  end
  effinal(k,1)=efxk;
  effinal(k,2)=efyk;
  end
  effinal