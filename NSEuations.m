%Solve Navier Stokes equation for simple shear flow
%SOR for Pressure boundary condition
clc
clear
Re = 10;                %Reynolds Number
%ft = 4 ;               %final time
lx = 1;                 %lenght of the domain
ly = 1;                 %width of the domain
nx = 5;                 %gridpoints along x
ny = 5;                 %gridpoints along y
dx = lx/(nx-1); dy = ly/(ny-1); %length of the grid cell

maxError = 0.001; maxit=200; nstep= 100; beta=1.2; dt=0.00125;
gx =0; gy=-100;                                   %Numerical variables

%---------------------------------------------
%Initial conditions
u = zeros(nx, ny-1); ut = zeros(nx,ny-1);  %velocity in x direction
v = zeros(nx-1, ny); vt = zeros(nx-1,ny);  %velocity in y direction
p = zeros(nx-1, ny-1);
tmp1 = zeros(nx-1, ny-1); tmp2 = zeros(nx-1, ny-1); %pressure
%set grid points
for i = 1:nx
    for j = 1:ny
          x=dx*(i-1); 
          y=dy*(j-0.5);
    end
end

%Initial velocities
uT = 1; vT = 0;            %at y=1
uB = -1; vB = 0;             %at y=0
%---------time loop----------%

%boundary conditions
u(1:nx,1) = 2*uB-u(1:nx,ny-2);
v(1,1:ny) = vT;
u(1:nx,ny-1) = 2*uT-u(1:nx,2);
v(nx-1,1:ny) = vT;
%--------------------------------------------------------------------
%--------------------------------------------------------------------
for is=1:nstep
% u-momentum
for i=2:nx-1
    for j = 2:ny-2
        F(i,j)= u(i,j)+dt*(-0.25*(((u(i+1,j)+u(i,j))^2-(u(i,j)+ ...
                u(i-1,j))^2)/dx+((u(i,j+1)+u(i,j))*(v(i,j)+ ...
                v(i-1,j))-(u(i,j)+u(i,j-1))*(v(i,j-1)+v(i-1,j-1)))/dy)+ ...
                (u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+ ...
                (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2+gx);
    end
end
% v-momentum 
for i=2:nx-2
    for j=2:ny-2
        G(i,j)=v(i,j)+dt*(-0.25*(((u(i,j+1)+u(i,j))*(v(i+1,j)+ ...
                v(i,j))-(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)))/dx+ ...
                ((v(i,j+1)+v(i,j))^2-(v(i,j)+v(i,j-1))^2)/dy)+ ...
                (v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+ ...
                (v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2 +gy);
    end
end
 %Solve poisson equation
 for i=2:nx-2
   for  j=2:ny-2
     p(i,j)= (1/dt)*((F(i,j)-F(i-1,j))/dx+((G(i,j)-G(i,j-1))/dy));
 end
 end
end
F
G

%Velocity update 


