function [w,cor_A,cor_phi]=allphase_fast(x,N,win,draw)
  
   if nargin< 4
draw=false;
  end

  
x=x(:);
%h=blackman(N);
ind=2*pi*[0:N-1]'/N;
%h=0.42-0.5*cos(ind)+0.08*cos(2*ind);
%X1=fft(x(1:N).*h);
%X2=fft(x(2:N+1).*h);


a=[   419.34299 ,...
  -407.33132,...
   281.22500,...
   -92.66868,...
     9.10350]/1000;
%a=[0.42,-0.5,0.08];

a=[   376.18589,...
  -369.22436,...
   287.02564,...
  -130.77564,...
    24.88141]/1000;
a=[0.42,-0.25,0.04];
a=win;
X1=fft(x(1:N));
X2=x(N+1)-x(1)+X1;
X2=X2.*exp(2*pi*1i*[0:N-1]'/N);

Z1=a(1)*X1;
Z2=a(1)*X2;
h=a(1);
for k=2:length(a)
  Z1=Z1+a(k)*circshift(X1,k-1)+a(k)*circshift(X1,-k+1);
  Z2=Z2+a(k)*circshift(X2,k-1)+a(k)*circshift(X2,-k+1);
  h=h+a(k)*cos(ind*(k-1))*2;
end
X1=Z1;
X2=Z2;


X=abs(X1);
if draw
figure;plot(X);
end
kstar=find( (X>[0;X(1:end-1)]) & (X>[0;0;X(1:end-2)]) &...
            (X>[X(2:end);0]) & (X>[X(3:end);0;0]) &  ...
            X>mean(X)*3);
%kstar=kstar(1:ceil(length(kstar)/2));
w=N/2/pi*angle(X2(kstar)./X1(kstar));
w=w+N*(1-sign(w))/2;

L=length(kstar);
phi=zeros(L,1);
A=zeros(L,1);
for k=1:L
  
  cor_w=w(k);
  
  y= exp(2*pi*1i*cor_w/N*[0:N-1]');
  Y=fft(y.*h);
  phi(k) = angle(X1(kstar(k))./Y(kstar(k)));
  newy=exp(2*pi*1i*cor_w/N*[0:N-1]'+phi(k)*1i);
  newY=fft(newy.*h);
  A(k) = abs(X1(kstar(k))/newY(kstar(k)));
end

cor_phi=phi;
cor_A=A;