%Stokes equation for 2D channel for 
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
MaxErr = 0.001;         %maximum error for pressure
Maxit=10;              %maximum iteration
nstep= 10;
beta=1.2;               %SOR factor
dt=0.01;               %time step size
gx =0; gy=0;        %external forces

%---------------------------------------------
%Variable initialization
u = zeros(nx+1,ny+2); ut = zeros(nx+1,ny+2);  %velocity in x direction
v = zeros(nx+2,ny+1); vt = zeros(nx+2,ny+1);  %velocity in y direction
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
%Check to ensure that the defined time step meets stability condition?
if dt>0.25*(min(dx,dy))^2/visc || dt>(1/Re)
    dt=min((0.25*(min(dx,dy))^2/visc),(visc));
end
%-----------------------------------------------------------------------
%Set boundary conditions
for n=1:nstep
    
    for j=2:ny+1
        y=dy*(j-1.5);
        u(1,j)= y*(1-y);               %u at left edge
        u(nx+1,j)=y*(1-y);            %u at right edge
        
    end
    
    for i=2:nx+1
        v(i,1)=0;                   %v at bottom edge
        v(i,ny+1)=0;                %v at top edge
    end
    
    
    for i=2:nx+1
        u(i,ny+2)=2-u(i,ny+1);       %u at top edge
        u(i,1)=-2-u(i,2);            %u at bottom edge
        
    end
    
    for j=2:ny
        v(nx+2,j)=-v(nx+1,j);        %v at right
        v(1,j)=-v(2,j);              %v at left
        
    end
    
    
%     for i=2:nx+1
%         p(i,ny+2)=6-p(i,ny+1);              %at top
%         p(i,1)=6-p(i,2);                    %at bottom
%     end
%     
%     for j=2:ny+1
%         p(nx+2,j)=6-p(nx+1,j);              %at right
%         p(1,j)=6-p(2,j);                    %at left
%     end
    
    %p
    %------------------------------------------------------------
    %temporary velocity
    %ut=u;
    for i=2:nx
        for j=2:ny+1
            ut(i,j)=u(i,j)+dt*(visc*((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+...
                (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2)+gx);
        end
    end
    
    %vt=v;
    for i=2:nx+1
        for j=2:ny
            vt(i,j)=v(i,j)+dt*(visc*((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+...
                (v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2)+gy);
        end
    end
    
    %ut
    %vt
    diag = (0.5/(1/dx^2+1/dy^2));
    for it=1:Maxit
        p_chk=p;
        for i=2:nx+1
            for j=2:ny+1
                %p(i,j)=3;
%                 rhs= (1/dt)*((ut(i,j)-ut(i-1,j))/dx...
%                     +(vt(i,j)-vt(i,j-1))/dy);
%                 p(i,j)= (1-beta)*p(i,j)+...
%                     beta*diag*( (1/dx^2)*(p(i+1,j)+p(i-1,j))+...
%                     (1/dy^2)*(p(i,j+1)+p(i,j-1))-rhs);
%                
                p(i,j)=  p(i,j) + ...
                          diag*( (1/dx^2)*(p(i+1,j)-2*p(i,j)+p(i-1,j))+ ...
                                 (1/dy^2)*(p(i,j+1)-2*p(i,j)+p(i,j-1)) ...
                -(Re/dt)*((ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy) );
            end
        end
        for i=2:nx+1
            p(i,ny+2)= p(i,ny+1);              %at top
            p(i,1)= p(i,2);                    %at bottom
        end     
        for j=2:ny+1
            p(nx+2,j)=2*(1-2*visc)-p(nx+1,j);   %at right
            p(1,j)=2*p(2,j);                    %at left
        end
        
       if max(max(abs(p_chk-p))) <MaxErr, break, end

    end
  
    
    
    %update velocities
    %u velocity
    
    for i=2:nx
        for j= 2:ny+1
            u(i,j)= ut(i,j)- ((dt/Re)*((p(i+1,j)-p(i,j))/dx));
        end
    end
    %v-velocity
    for i=2:nx+1
        for j= 2:ny
            v(i,j)= vt(i,j)- ((dt/Re)*((p(i,j+1)-p(i,j))/dy));
        end
    end
    u
    v
    p
   
end






        
        


        