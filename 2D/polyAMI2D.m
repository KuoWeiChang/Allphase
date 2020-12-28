function [efx,efy]=polyAMI2D(x,N,nx,ny,fx,fy,dirc)
  if dirc==1
    xe=exp_maker2(fx,fy,1,0,N,nx,ny);
    xb=exp_maker2(fx-0.5,fy,1,0,N,nx,ny);
    xf=exp_maker2(fx+0.5,fy,1,0,N,nx,ny);
     g=[fx-0.5;fx;fx+0.5];
  else
     xe=exp_maker2(fx,fy,1,0,N,nx,ny);
    xb=exp_maker2(fx,fy-0.5,1,0,N,nx,ny);
    xf=exp_maker2(fx,fy+0.5,1,0,N,nx,ny);
     g=[fy-0.5;fy;fy+0.5];
  end
b=abs([sum(sum(conj(xb).*x));sum(sum(conj(xe).*x));sum(sum(conj(xf).*x))]);
if dirc==1
  efx=-1/2*(g(1)^2*(b(2)-b(3))+g(2)^2*(b(3)-b(1))+g(3)^2*(b(1)-b(2)))/(g(1)*(b(3)-b(2))+g(2)*(b(1)-b(3))+g(3)*(b(2)-b(1)));
  efy=fy;
else
  efy=-1/2*(g(1)^2*(b(2)-b(3))+g(2)^2*(b(3)-b(1))+g(3)^2*(b(1)-b(2)))/(g(1)*(b(3)-b(2))+g(2)*(b(1)-b(3))+g(3)*(b(2)-b(1)));;
  efx=fx;
end
  