function s=exp_maker2(fx,fy,A,phi,N,nx,ny)
x=0;
W=exp(-2*pi*1i/N);
for k=1:length(fx)
  x=x+A(k)*exp(2*pi*1i*ny/N*fy(k))*exp(2*pi*1i*nx/N*fx(k))*exp(1i*phi(k));
end
s=x;

