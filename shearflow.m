%Divergent test for simple shear flow and channel flow
clc
clear
Re = 10;       %Reynolds Number
dt = .01;       %Time step
%ft = 4 ;        %final time
lx = 1;         %lenght of the domain
ly = 1;         %width of the domain
nx = 5;        %gridpoints along x
ny = 5;        %gridpoints along y
dx = lx/(nx-1); dy = ly/(ny-1); %length of the grid cell
%---------------------------------------------
%Initial conditions
u = zeros(nx, ny-1); um = zeros(nx,ny-1);  %velocity in x direction
v = zeros(nx-1, ny); vm = zeros(nx-1,ny);  %velocity in y direction
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
          %u(i,j)=-y*(y-1);
    end
end

for i = 1:nx-1
    for j = 1:ny
          x=dx*(i-1);       
          y=dy*(j-0.5);    
          v(i,j) = 0;
    end
end

%ub = -y(j)*(y(j)-1), vb =0;
for i=1:nx-1
    for j =1:ny-1
        d(i,j)= (u(i+1,j)-u(i,j))/dx + (v(i,j+1)-v(i,j))/dy;
    end
end

quiver(u,v)

%--------------------------------------------------------------------
%--------------------------------------------------------------------
% u-momentum
for i=2:nx-1
    for j = 2:ny-2
        um(i,j)= (p(i,j)-p(i-1,j))/dx + (1/Re)*(((u(i+1,j)+u(i-1,j)...
        -2*u(i,j))/dx^2 + (u(i,j+1)+u(i,j-1)-2*u(i,j))/dy^2));
    end
end
% v-momentum 
for i=2:nx-2
    for j=2:ny-1
        vm(i,j)= (p(i,j)-p(i,j-1))/dy + (1/Re)*(((v(i+1,j)+v(i-1,j)...
        -2*v(i,j))/dx^2 + (v(i,j+1)+v(i,j-1)-2*v(i,j))/dy^2));
    end
end



        




