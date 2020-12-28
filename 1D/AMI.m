function [f,amp,ph]=AMI(x,iter,pm,ef)
  x=x(:);
  N=length(x);
  if nargin< 4
    X=fft(x);
    [M,k0]=max(abs(X));
    ef=k0-1;
  end

  for m=1:iter
    s1=exp(-2*pi*1i*[0:N-1]*(ef+pm)/N);
    s2=exp(-2*pi*1i*[0:N-1]*(ef-pm)/N);
    S1=abs(s1*x);
    S2=abs(s2*x);
    h=1/2*(S1-S2)/(S1+S2);
    ef=ef+h;
  end
  h=blackman(N);
  
  y= exp(2*pi*1i*ef/N*[0:N-1]');
  p=x(:);
  Y=fft(y.*h);
  X=fft(p.*h);
  n=round(ef)+1;
  phi= angle(X(n)./Y(n));
  newy=exp(2*pi*1i*ef/N*[0:N-1]'+phi*1i);
  newY=fft(newy.*h);
  A= abs(X(n)./newY(n));
  
 f=ef;
 ph=phi;
 amp=A;