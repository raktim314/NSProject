%unsteady stokes equation for Lid-Driven Cavity flow
%SOR method for faster convergent
clc;
clear all;
close all;
%-----------------------------------------------------------------------
%Define all numerial values
%------------------------------------------------------------------------
Re = 1;                 %Reynolds Number
visc=1/Re;              %viscosity
nx = 64;                 %gridpoints along x
ny = 64;                 %gridpoints along y
lx = 1;                 %lenght of the domain
ly = 1;                 %width of the domain
ft=0.5;                 %final time
MaxErr = 0.001;       %maximum error for pressure
Maxit=100;              %maximum iteration
nstep= 1000;
beta=1.2;               %SOR factor
dt=0.01;               %time step size
gx =0; gy=0;        %external forces

%---------------------------------------------
%Variable initialization
u = zeros(nx+1,ny+2); ut = zeros(nx+1,ny+2);  %velocity in x direction
v = zeros(nx+2,ny+1); vt = zeros(nx+2,ny+1);  %velocity in y direction
p = zeros(nx+2,ny+2);                       %pressure
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
   
    for j=2:ny+1
        u(1,j)=0;               %u at left edge
        u(nx+1,j)=0;            %u at right edge
    end
    
     for i=2:nx+1
        u(i,ny+2)=2-u(i,ny+1);       %u at top edge
        u(i,1)= -u(i,2);             %u at bottom edge
    end
    
    for i=2:nx+1
        v(i,1)=0;                   %v at bottom edge
        v(i,ny+1)=0;                %v at top edge
    end
    
    for j=2:ny
        v(nx+2,j)=-v(nx+1,j);        %v at right
        v(1,j)=-v(2,j);              %v at left
    end
    
    
    for i=2:nx+1
        p(i,ny+2)= p(i,ny+1);              %at top
        p(i,1)= -p(i,2);                    %at bottom
    end
    for j=2:ny+1
        p(nx+2,j)= p(nx+1,j);              %at right
        p(1,j)= p(2,j);                    %at left
    end
    
   %------------------------------------------------------------
    for n=1:nstep
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
    %-------------------------------------------------------------------
    %boundary velocities update
    
    for j=2:ny+1
        ut(1,j)=0;               %u at left edge
        ut(nx+1,j)=0;            %u at right edge
    end
    
    for i=2:nx+1
        ut(i,ny+2)=2-ut(i,ny+1);       %u at top edge
        ut(i,1)= -ut(i,2);             %u at bottom edge
    end
    
    for i=2:nx+1
        vt(i,1)=0;                   %v at bottom edge
        vt(i,ny+1)=0;                %v at top edge
    end
    
    for j=2:ny
        vt(nx+2,j)=-vt(nx+1,j);        %v at right
        vt(1,j)=-vt(2,j);              %v at left
    end
    %ut
    %vt
     for i=2:nx+1
        p(i,ny+2)= p(i,ny+1);              %at top
        p(i,1)= -p(i,2);                    %at bottom
    end
    for j=2:ny+1
        p(nx+2,j)= p(nx+1,j);              %at right
        p(1,j)= p(2,j);                    %at left
    end
    %------------------------------------------------------------------
    %Solve Pressure
    %------------------------------------------------------------------
      diag = (0.5/(1/dx^2+1/dy^2));
%     rhs= (1/dt)*((ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy);
    for it=1:Maxit
        p_chk=p;
        for i=2:nx+1
            for j=2:ny+1
%                 p(i,j)= (1-beta)*p(i,j)+...
%                     beta*diag*( (1/dx^2)*(p(i+1,j)+p(i-1,j))+...
%                     (1/dy^2)*(p(i,j+1)+p(i,j-1))-rhs);
               
                p(i,j)= p(i,j) + ...
                          diag*( (1/dx^2)*(p(i+1,j)-2*p(i,j)+p(i-1,j))+ ...
                                 (1/dy^2)*(p(i,j+1)-2*p(i,j)+p(i,j-1)) ...
                -(Re/dt)*((ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy) );
            end
        end
        
        Err=0.0;
        for i=2:nx+1
            for j=2:ny+1
                Err=Err+abs(p_chk(i,j)-p(i,j));
            end
        end
        if Err <= MaxErr, break, end
        
        for i=2:nx+1
            p(i,ny+2)= p(i,ny+1);              %at top
            p(i,1)= -p(i,2);                    %at bottom
        end     
        for j=2:ny+1
            p(nx+2,j)= p(nx+1,j);              %at right
            p(1,j)= p(2,j);                    %at left
        end
      

    end
%------------------------------------------------------------------------ 
%update velocities
%------------------------------------------------------------------------
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
 
    end
t =t+dt
%relocate the grid points
u_avr(1:nx+1,1:ny+1)= 1/2 *(u(1:nx+1,1:ny+1)+u(1:nx+1,2:ny+2));
v_avr(1:nx+1,1:ny+1)= 1/2 *(v(1:nx+1,1:ny+1)+v(2:nx+2,1:ny+1));
p_avr(1:nx+1,1:ny+1)=1/4*(p(1:nx+1,1:ny+1)+p(2:nx+2,1:ny+1)...
    +p(1:nx+1,2:ny+2)+p(2:nx+2,2:ny+2));
wt(1:nx+1,1:ny+1)=(v(2:nx+2,1:ny+1)-v(1:nx+1,1:ny+1))/dx...
    -(u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1))/dy;
x(1:nx+1)=(0:nx);
y(1:ny+1)=(0:ny);
% %-------------------------------------------------------------------------%
% %Plot the variables
% %-------------------------------------------------------------------------%
% figure(1), quiver(x,y,(rot90(fliplr(u_avr))),(rot90(fliplr(v_avr)))),...
% xlabel('nx'),ylabel('ny'),title('Velocity Vectour Plot for Grid size 16x16');
% axis([0 nx 0 ny]),axis
% ('square');
% figure(2),
% subplot(211), contourf(x,y,rot90(fliplr(u_avr)),30),...
% xlabel('nx'),ylabel('ny'),title('u-velocity for 32x32 grid'); axis('square','tight'); colorbar
% subplot(212), contourf(x,y,rot90(fliplr(v_avr)),30),...
% xlabel('nx'),ylabel('ny'),title('v-velocity for 32x32 grid'); axis('square','tight'); colorbar
% figure(3),
% subplot(211), contourf(x,y,rot90(fliplr(p_avr)),30),...
% xlabel('nx'),ylabel('ny'),title('Pressure for 32x32 grid'); axis('square','tight'); colorbar
% subplot(212), contourf(x,y,rot90(fliplr(wt)),30),...
% xlabel('nx'),ylabel('ny'),title('Vorticity for 32x32 grid') ;axis('square','tight'); colorbar
% figure(1), quiver(x,y,(rot90(fliplr(u_avr))),(rot90(fliplr(v_avr)))),...
% xlabel('nx'),ylabel('ny'),title('Velocity Vectour Plot');
%axis([0 nx 0 ny]),axis('square');
 figure(2), contour(x,y,(rot90(fliplr(p_avr)))),...
     title('pressue for grid size 32x32 '); colorbar
% figure(3), plot(u(floor(nx/2),2:ny+1),2:ny+1, 'or-'),title('u velocity');
% xlabel('u')
% ylabel('y')
% hold on
% figure(4), plot(1:nx+1,v(1:nx+1, floor(ny/2)), 'ob-'),title('v-velocity')
% xlabel('x')
% ylabel('v')
% hold on;
        


        