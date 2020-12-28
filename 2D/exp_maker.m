function s=exp_maker(f,A,phi,N,n)
x=0;
W=exp(-2*pi*1i/N);
for k=1:length(f)
  x=x+A(k)*exp(2*pi*1i*n/N*f(k)+1i*phi(k));
end
s=x;

