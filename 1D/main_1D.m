%% test allphase fast
clear all;
close all;
clc;

N=400;

%exp1 1D all phase

x=exp_maker([128.487,23.31],[1.2,0.8],[0.311,0.911],N,[0:N]');
%x=x+randn(N+1,1)*0.3;
%x=x+randn(length(x),1)*0.5;
[Fe,Ae,phie]=allphase_fast(x,N,[0.42,-0.25,0.04])

[Estf,EA,Ephi]=AMI(x(1:N),5,0.5,128)
[Estf2,EA2,Ephi2]=AMI(x(1:N),5,0.5,23)
 
 