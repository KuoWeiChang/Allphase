%% 1 line vertical and 2 lines horizon
function [fx,fy,A,phi]=AMI_2D_all(x,N,nx,ny,method,h)
  v1=x(:,1);
  h1=x(1,:);
  h2=x(2,:);
  
  [ef1,eA1,ephi1]=AMI_all(h1,N,nx,method,h);
    [ef2,eA2,ephi2]=AMI_all(h2,N,nx,method,h);
    vf=(ephi2-ephi1)*N/2/pi;
    L=length(vf);
    
        [vfr,vAr,vphir]=AMI_all(v1,N,ny,method,h);

    evf=zeros(L,1);
   for k=1:L
     [tmp,ind]=min(abs(vf(k)-vfr));  %% find nearest neighbor and substitute
      evf(k)=vfr(ind);
   end
   
    fx=ef1;
    fy=evf;
    A=eA1;
    phi=ephi1;