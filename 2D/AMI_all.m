function [f,A,phi]=AMI_all(x,N,n,method,h)
  
    if nargin<5
    h=ones(N,1);
  end
  x=x(:);
  h=h(:);
  n=n(:);
  
  X=abs(fft(x));
  %% step 1 : find peaks
  ind = find( (X>[0;X(1:end-1)]) & (X>[0;0;X(1:end-2)]) &...
            (X>[X(2:end);0]) & (X>[X(3:end);0;0]) &  ...
            X>mean(X)*3);

  xin = x.*h;
  z=zeros(length(ind),3);
  B=zeros(N,length(ind));
  for k=1:length(ind)
    ef=ind(k)-1;
    if method == 'polyAMI2'
      for iter=1:5
        ef=polyAMI2(xin,N,n,ef-0.5,ef,ef+0.5);
      end
    end
    z(k,1)=ef;
    %{
    y=exp(2*pi*1i*n*ef/N);
    a=y'*xin;
    b=y'*(y.*h);
    pha = angle(a/b);  
    amp = abs(a/b);
    z(k,2)=amp;
    z(k,3)=pha;
  %}
  B(:,k)=exp_maker(ef,1,0,N,n);
  
  end

r=pinv(B)*x;


  f=z(:,1);
  A=abs(r);
  phi=angle(r);