%Solve NS equation for simple shear flow
clc
clear
Re = 10;       %Reynolds Number
dt = .001;       %Time step
%ft = 4 ;        %final time
lx = 1;         %lenght of the domain
ly = 1;         %width of the domain
nx = 5;        %gridpoints along x
ny = 5;        %gridpoints along y
dx = lx/(nx-1); dy = ly/(ny-1); %length of the grid cell
%---------------------------------------------
%Initial conditions
u = zeros(nx, ny-1); ut = zeros(nx,ny-1);  %velocity in x direction
v = zeros(nx-1, ny); vt = zeros(nx-1,ny);  %velocity in y direction
p = zeros(nx-1, ny-1);
px = zeros(nx-1, ny-1); py = zeros(nx-1, ny-1); %pressure
%set grid points
for i = 1:nx
    for j = 1:ny-1
          x=dx*(i-1); 
          y=dy*(j-0.5);
    end
end
[X,Y] = meshgrid(y,x);
% u-velocity 
for i = 1:nx
    for j = 1:ny-1
          x=dx*(i-1);
          y=dy*(j-0.5);
          u(i,j) = 2*y-1;
       
    end
end

for i = 1:nx-1
    for j = 1:ny
          x=dx*(i-1);       
          y=dy*(j-0.5);    
          v(i,j) = 0;
    end
end

%----------------------------------------------------
%Initial velocities
uT = 1; vT = 0;            %at y=1
uB = -1; vB = 0;             %at y=0
%---------time loop----------%

for i = 1:100
%boundary conditions
u(1:nx,ny-1) = 2*uB-u(1:nx,ny-2);
v(1:nx-1,ny) = vT;
u(1:nx,1) = 2*uT-u(1:nx,2);
v(nx-1,1:ny-1) = vT;

for i = 2:nx-1
    for j = 2:ny-2
        ut(i,j)= u(i,j)+dt*((1/Re)*(((u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2)...
        +((u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2))-((u(i+1,j)+u(i,j))^2     ...
        -(u(i,j)+u(i-1,j))^2)/dx-((u(i,j+1)+u(i,j))*(v(i,j)+v(i-1,j))   ...
        -(u(i,j)+u(i,j-1))*(v(i,j-1)+v(i-1,j-1)))/dy); %p constant and dp/dx=0
    end
end

for i=2:nx-2
    for j=2:ny-1
        vt = v(i,j)+dt*((1/Re)*(((v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2)...
        +((v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2))-(((u(i,j)+u(i,j-1))*(...
        v(i+1,j)+v(i,j))-(u(i-1,j)+u(i-1,j-1))*v(i,j)+v(i-1,j)))/dx...
        -((v(i,j+1)+v(i,j))^2-(v(i,j)+v(i,j-1))^2)/dy);%p is constant &dp/dy=0
    end
end
end
ut
vt

