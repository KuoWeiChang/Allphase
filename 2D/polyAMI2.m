function f=polyAMI2(x,N,n,fb,fe,ff)
  xe=exp(2*pi*1i*n*fe/N);
  xb=exp(2*pi*1i*n*fb/N);
  xf=exp(2*pi*1i*n*ff/N);
  g=[fb;fe;ff];
 % A=[g.^2,g,ones(3,1)];
 % det(A)
  b=[abs(xb'*x);abs(xe'*x);abs(xf'*x)];
  %v=inv(A)*b;
  %f=v(2)/v(1)/-2;
  f=-1/2*(g(1)^2*(b(2)-b(3))+g(2)^2*(b(3)-b(1))+g(3)^2*(b(1)-b(2)))/(g(1)*(b(3)-b(2))+g(2)*(b(1)-b(3))+g(3)*(b(2)-b(1)));
