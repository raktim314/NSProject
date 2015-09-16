
%unsteady Navier-stokes equation for simple shear flow
%Solve pressure with projection method
%SOR method for faster convergent
%External forces are considered to be zero
clc;
clear all;
close all;
%-----------------------------------------------------------------------
%Define all numerial values
%------------------------------------------------------------------------
Re = 1;                     %Reynolds Number
visc=1/Re;                  %viscosity
nx = 16;                    %gridpoints along x
ny = 16;                    %gridpoints along y
lx = 1;                     %lenght of the domain
ly = 1;                     %width of the domain
ft=0.5;                     %final time
MaxErr = 0.001;             %maximum error for pressure
Maxit=10;                   %maximum iteration
nstep= 100;                 %time step
beta=1.2;                   %SOR factor
dt=0.01;                    %time step size


%---------------------------------------------
%Variable initialization
u = zeros(nx+1,ny+2); ut = zeros(nx+1,ny+2);    %velocity in x direction
v = zeros(nx+2,ny+1); vt = zeros(nx+2,ny+1);    %velocity in y direction
p = zeros(nx+2,ny+2);                           %pressure
t=0;                                            %Initial time
%----------------------------------------------------------------------
%set grid points
dx=lx/nx; dy=ly/ny;                             %length of the grid cell
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
%Time integration loop start
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
        x=dx*(i-1.5);
        v(i,1)=0;
        v(i,ny+1)=0;
    end
    
%-------------------------------------------------------------------------
    %Set boundary conditions
    for i=2:nx+1
        u(i,ny+2)=2-u(i,ny+1);              %at top line
        u(i,1)=-2-u(i,2);                   %at bottom line
    end
    
    for j=2:ny
        v(nx+2,j)=-v(nx+1,j);               %at right
        v(1,j)=-v(2,j);                     %at left
    end
    
    
    for i=2:nx+1
        p(i,ny+2)=6-p(i,ny+1);              %at top
        p(i,1)=6-p(i,2);                    %at bottom
    end
    
    for j=2:ny+1
        p(nx+2,j)=6-p(nx+1,j);              %at right
        p(1,j)=6-p(2,j);                    %at left
    end
    
    %p
    %------------------------------------------------------------
    %temporary velocity
    %ut=u;
    for i=2:nx
        for j=2:ny+1
            usqr_1=(1/2*(u(i+1,j)+u(i,j)))^2;
            usqr_2=(1/2*(u(i-1,j)+u(i,j)))^2;
            uv_x1=1/2*(u(i,j)+u(i,j+1))*1/2*(v(i,j)+v(i+1,j));
            uv_x2=1/2*(u(i,j)+u(i,j-1))*1/2*(v(i,j-1)+v(i+1,j-1));
            A=(usqr_1-usqr_2)/dx+(uv_x1-uv_x2)/dy; %Advection term X-direction
            D2_u=((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2)+...
                ((u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2); %Laplacian operator for u
            ut(i,j)=u(i,j)+dt*(-A+visc*D2_u);
          
        end
    end
    
    %vt=v;
    for i=2:nx+1
        for j=2:ny
            vsqr_1=(1/2*(v(i,j+1)+v(i,j)))^2;
            vsqr_2=(1/2*(v(i,j)+v(i,j-1)))^2;
            uv_y1=1/2*(u(i,j+1)+u(i,j))*1/2*(v(i,j)+v(i+1,j));
            uv_y2=1/2*(u(i-1,j+1)+u(i-1,j))*1/2*(v(i,j)+v(i-1,j));
            B=(uv_y1-uv_y2)/dx+(vsqr_1-vsqr_2)/dy; %Advection term Y-direction
            D2_v=((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2)+...
                ((v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2); %Laplacian operator for v
            vt(i,j)=v(i,j)+dt*(-B+visc*D2_v);
            
        end
    end
    
    diag = (0.5/(1/dx^2+1/dy^2));
    for it=1:Maxit
        p_chk=p;
        for i=2:nx+1
            for j=2:ny+1
                p(i,j)= p(i,j) + ...
                          diag*( (1/dx^2)*(p(i+1,j)-2*p(i,j)+p(i-1,j))+ ...
                                 (1/dy^2)*(p(i,j+1)-2*p(i,j)+p(i,j-1)) ...
                -(Re/dt)*((ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy) );
               
            end
        end
        for i=2:nx+1
            p(i,ny+2)=6-p(i,ny+1);              %at top
            p(i,1)=6-p(i,2);                    %at bottom
        end     
        for j=2:ny+1
            p(nx+2,j)=6-p(nx+1,j);              %at right
            p(1,j)=6-p(2,j);                    %at left
        end
%check pressure convergence 
        Err=0.0;
        for i=2:nx+1
            for j=2:ny+1
                Err=Err+abs(p_chk(i,j)-p(i,j));
            end
        end
if Err <= MaxErr, break, end
    end
%-----------------------------------------------------------------------  
%Update velocities
%Find the final velocity
%-----------------------------------------------------------------------
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
end                             %Time integration end 






  

        

























