%unsteady stokes equation for simple shear flow 
%pressure is not constant
%SOR method for faster convergent 
clc;
clear all;
close all;
%-----------------------------------------------------------------------
%Define all numerial values
%------------------------------------------------------------------------
Re = 1;                 %Reynolds Number
visc=1/Re;              %viscosity
nx = 4;                 %gridpoints along x
ny = 4;                 %gridpoints along y
lx = 1;                 %lenght of the domain
ly = 1;                 %width of the domain
ft=0.5;                 %final time
maxError = 0.001;       %maximum error for pressure
maxit=100;              %maximum iteration
nstep= 1; 
beta=1.2;               %SOR factor
dt=0.001;               %time step size
gx =0; gy=0;        %external forces                           

%---------------------------------------------
%Variable initialization
u = zeros(nx+1,ny+2); ut = zeros(nx,ny-1);  %velocity in x direction
v = zeros(nx+2,ny+1); vt = zeros(nx-1,ny);  %velocity in y direction
p = zeros(nx+2,ny+2);                       %pressure
tmp1 = zeros(nx+2,ny+2); tmp2 = zeros(nx+2,ny+2); %pressure
t=0;                                                %Initial time
%----------------------------------------------------------------------
%set grid points
dx=lx/nx; dy=ly/ny;                         %length of the grid cell
for i = 0:nx+2
    for j = 0:ny+2
          x=dx*(i-1); 
          y=dy*(j-1.5);
    end
end
%-----------------------------------------------------------------------
%Interior nodes
for n=1:nstep
%set u
for j=2:ny+1
    y=dy*(j-1.5);
    u(1,j)=2*y-1;
    u(nx+1,j)=2*y-1;
end
%set v
for i=2:nx+1
    x=dx*(i-1);
    v(i,1)=0;
    v(i,ny+1)=0;
end
%set pressure p (assume pressure is constant & =1)
for i=2:nx+1
    for j=2:ny+1
        p(i,j)=1;
    end
end
%-------------------------------------------------------------------------
%Set boundary conditions
for i=2:nx+1
    u(i,ny+1)=2-u(i,nx+2);              %at top line
    u(i,1)=-2-u(i,2);                   %at bottom line
end

for j=2:ny
    v(nx+2,j)=-v(nx+1,j);               %at right
    v(1,j)=-v(2,j);                     %at left
end

for i=2:nx+1
   for j=2:ny+1
    p(nx+2,j)=6-p(nx+1,j);              %at right
    p(i,ny+2)=6-p(i,ny+1);              %at top
    p(1,j)=6-p(2,j);                    %at left
    p(i,1)=6-p(i,2);                    %at bottom
   end
end
%------------------------------------------------------------
%temporary velocity
for i=2:nx
    for j=2:ny+1
        ut(i,j)=u(i,j)+dt*(visc*((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+...
            (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2)+gx);
    end
end

for i=2:nx+1
    for j=ny
        vt(i,j)=v(i,j)+dt*(visc*((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+...
            (v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2)+gy);
    end
end

for it=1:maxit
    p_chk=p;
    for i=2:nx+1
        for j=2:ny+1
            p(i,j)=(beta/2)*( (1/(1/dx^2+1/dy^2))*((1/dx^2)*(p(i+1,j)...
                +p(i-1,j))+ (1/dy^2)*(p(i,j+1)+p(i,j-1))-...
                (Re/dt)*((ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy)))+...
                (1-beta)*p(i,j);
        end
    end
    if max(max(abs(p_chk-p))) <maxError, break, end
end
end

%update velocities
%u velocity 

for i=2:nx-1
    for j= 2:ny-2
        u(i,j)= ut(i,j)- ((dt/Re)*((p(i+1,j)-p(i,j))/dx));
    end
end
%v-velocity 
for i=2:nx-2
    for j= 2:ny-1
        v(i,j)= vt(i,j)- ((dt/Re)*((p(i,j+1)-p(i,j))/dy));
    end
end





        
        


        