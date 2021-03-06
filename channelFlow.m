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
nx = 64;                 %gridpoints along x
ny = 64;                 %gridpoints along y
lx = 1;                 %lenght of the domain
ly = 1;                 %width of the domain
ft=0.5;                 %final time
MaxErr = 0.001;         %maximum error for pressure
Maxit=10;              %maximum iteration
nstep= 800;
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
    y=dy*(j-1.5);
    u(1,j)= y*(1-y);               %u at left edge
    u(nx+1,j)=y*(1-y);            %u at right edge
    
end

for i=2:nx+1
    v(i,1)=0;                   %v at bottom edge
    v(i,ny+1)=0;                %v at top edge
end


for i=2:nx+1
    u(i,ny+2)= -u(i,ny+1);       %u at top edge
    u(i,1)= -u(i,2);            %u at bottom edge
    
end

for j=2:ny
    v(nx+2,j)=-v(nx+1,j);        %v at right
    v(1,j)=-v(2,j);              %v at left
    
end
for i=2:nx+1
    p(i,ny+2)= p(i,ny+1);              %at top
    p(i,1)= p(i,2);                    %at bottom
end
for j=2:ny+1
    p(nx+2,j)=2*(1-2*visc)-p(nx+1,j);   %at right
    p(1,j)=2-p(2,j);                    %at left
end

%---------------------------------------------------------------------
%Temporary velocity
%---------------------------------------------------------------------
for n=1:nstep
    for i=2:nx
        for j=2:ny+1
            ut(i,j)=u(i,j)+dt*(visc*((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+...
                (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2)+gx);
        end
    end
    
    
    for i=2:nx+1
        for j=2:ny
            vt(i,j)=v(i,j)+dt*(visc*((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+...
                (v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2)+gy);
        end
    end
    %---------------------------------------------------------------------
    %Update boundary values
    %---------------------------------------------------------------------
    for j=2:ny+1
        y=dy*(j-1.5);
        ut(1,j)= y*(1-y);               %u at left edge
        ut(nx+1,j)=y*(1-y);            %u at right edge
    end
    
    for i=2:nx+1
        vt(i,1)=0;                   %v at bottom edge
        vt(i,ny+1)=0;                %v at top edge
    end
    
    for i=2:nx+1
        ut(i,ny+2)= -ut(i,ny+1);       %u at top edge
        ut(i,1)= -ut(i,2);            %u at bottom edge
    end
    
    for j=2:ny
        vt(nx+2,j)=-vt(nx+1,j);        %v at right
        vt(1,j)=-vt(2,j);              %v at left
    end
    for i=2:nx+1
        p(i,ny+2)= p(i,ny+1);              %at top
        p(i,1)= p(i,2);                    %at bottom
    end
    for j=2:ny+1
        p(nx+2,j)=2*(1-2*visc)-p(nx+1,j);   %at right
        p(1,j)=2-p(2,j);                    %at left
    end
    %ut
    %vt
    %-----------------------------------------------------------------------
    %Solving Pressure
    %-----------------------------------------------------------------------
    diag = (0.5/(1/dx^2+1/dy^2));
    for it=1:Maxit
        p_chk=p;
        for i=2:nx+1
            for j=2:ny+1
                %                 rhs= (1/dt)*((ut(i,j)-ut(i-1,j))/dx...
                %                     +(vt(i,j)-vt(i,j-1))/dy);
                %                         p(i,j)= (1-beta)*p(i,j)+...
                %                             beta*diag*( ((1/dx^2)*(p(i+1,j)+p(i-1,j))+...
                %                             (1/dy^2)*(p(i,j+1)+p(i,j-1)))-rhs);
                p(i,j)=  p(i,j) + ...
                    diag*( (1/dx^2)*(p(i+1,j)-2*p(i,j)+p(i-1,j))+ ...
                    (1/dy^2)*(p(i,j+1)-2*p(i,j)+p(i,j-1)) ...
                    -(Re/dt)*((ut(i,j)-ut(i-1,j))/dx...
                    +(vt(i,j)-vt(i,j-1))/dy) );
            end
        end
        if max(max(abs(p_chk-p))) <MaxErr, break, end
        
        for i=2:nx+1
            p(i,ny+2)= p(i,ny+1);              %at top
            p(i,1)= p(i,2);                    %at bottom
        end
        for j=2:ny+1
            p(nx+2,j)=2*(1-2*visc)-p(nx+1,j);   %at right
            p(1,j)=2-p(2,j);                    %at left
        end
        
        
    end
    
    
    %----------------------------------------------------------------------
    %update velocity fields
    %----------------------------------------------------------------------
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
    %--------------------------------------------------------------------
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
        u(i,ny+2)= -u(i,ny+1);       %u at top edge
        u(i,1)= -u(i,2);            %u at bottom edge
    end
    
    for j=2:ny
        v(nx+2,j)=-v(nx+1,j);        %v at right
        v(1,j)=-v(2,j);              %v at left
        
    end
    %     u
    %     v
    %     p
end
%Error Calculation
e_u=0;
for i=2:nx
    for j=2:ny+1
        x=dx*(i-2);
        y=dy*(j-1.5);
        u_exact=y*(1-y);
        e_u=max(e_u, abs(u_exact-u(i,j)));
    end
end
e_v=0;
for i=2:nx+1
    for j=2:ny
        x=dx*(i-1.5);
        y=dy*(j-1);
        v_exact=0;
        e_v=max(e_v, abs(v_exact-v(i,j)));
    end
end
e_p=0;
for i=2:nx+1
    for j=2:ny+1
        x=dx*(i-1.5);
        y=dy*(j-1.5);
        p_exact=1-2*visc*x;
        e_p=max(e_p, abs(p_exact-p(i,j)));
    end
end
e_u
e_v
e_p

% %-----------------------------------------------------------------------
%relocate the grid points
%-----------------------------------------------------------------------
u_avr(1:nx+1,1:ny+1)= 0.5 *(u(1:nx+1,1:ny+1)+u(1:nx+1,2:ny+2));
v_avr(1:nx+1,1:ny+1)= 0.5 *(v(1:nx+1,1:ny+1)+v(2:nx+2,1:ny+1));
p_avr(1:nx+1,1:ny+1)=0.25*(p(1:nx+1,1:ny+1)+p(2:nx+2,1:ny+1)...
    +p(1:nx+1,2:ny+2)+p(2:nx+2,2:ny+2));
wt(1:nx+1,1:ny+1)=(v(2:nx+2,1:ny+1)-v(1:nx+1,1:ny+1))/dx...
    -(u(1:nx+1,2:ny+2)-u(1:nx+1,1:ny+1))/dy;
x(1:nx+1)=(0:nx);
y(1:ny+1)=(0:ny);
%-------------------------------------------------------------------------%
%Plot the variables
%-------------------------------------------------------------------------%
%figure(1), quiver(x,y,(rot90(fliplr(u_avr))),(rot90(fliplr(v_avr)))),...
%xlabel('nx'),ylabel('ny'),title('Velocity Vectour Plot');
%axis([0 nx 0 ny]),axis('square');
figure(2), contour(x,y,(rot90(fliplr(p_avr)))); colorbar
figure (3), plot(u(floor(nx/2),2:ny+1),2:ny+1,'p-');
xlabel('u-velocity');
ylabel('grid');
% hold on;
% figure(4), plot(u(floor(nx/2),2:ny+1),2:ny+1, '--');
hold on
hold all