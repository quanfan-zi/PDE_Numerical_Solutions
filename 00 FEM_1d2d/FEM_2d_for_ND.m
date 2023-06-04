function FEM_2d_for_ND
%% 
% u= cos(pi.*x).*cos(pi.*y);
% -\laplace u=2*pi*pi*cos(pi.*x).*cos(pi.*y);
% top and down  Dirichlet boundary condtion
% left and right Neumman boundary condtion
tic
format long;
pde.start_point_x=-1;
pde.end_point_x=1;
pde.start_point_y=-1;
pde.end_point_y=1;

N=5;
[L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(pde,N)

%M=20;
%[sDof,L2_error,H1_error]=run_main(M,pde)

toc
end

function [sDof,L2_error,H1_error]=run_main(M,pde)
%% intput M pde 维数
%% output sDof,L2_error,H1_error 维数

boundary_method=1;  % 1  2  3
domainLenght=pde.end_point_x-pde.start_point_x;
hmax=domainLenght/M;

% Generate mesh
[p,e,t]=meshGeneration(pde,M,hmax);
sDof=size(p,2);

% Generate stiffness matrix
A=stiffnessMatrix(p,t);

% Generate right hands
F=rightHands(A,p,e,t);

% boundary condition
 [A,F]=boundaryMatrix(A,F,p,e,boundary_method);

% solve A\F
u=solveAF(A,F,p,e,boundary_method);

% error estimate
[L2_error,H1_error]=errorEstimate(u,p,t);

% plot figure
plotFigure(p,e,t,u);




end

function [p,e,t]=meshGeneration(pde,M,hmax)
%% Generate mesh
g=[ 2     2     2     2
    -1     1     1    -1
     1     1    -1    -1
    -1    -1     1     1
    -1     1     1    -1
     1     1     1     1
     0     0     0     0];
[p,e,t]=initmesh(g,'Hmax',hmax);
%[p,e,t]=uniform_trimesh_delaunay(pde.start_point_x,pde.end_point_x,pde.start_point_y,pde.end_point_y,M,M);
end

function A=stiffnessMatrix(p,t)
%% Generate stiffness matrix
sDof=size(p,2);
A=sparse(sDof,sDof);
nel=size(t,2);
for iel=1:nel
    nd=t(1:3,iel);
    a=p(:,nd(1));
    b=p(:,nd(2));
    c=p(:,nd(3));
    k=elK(a,b,c);
    A(nd,nd)=A(nd,nd)+k;
end
end

function F=rightHands(A,p,e,t)
%% Generate right hands
sDof=size(p,2);
F=zeros(sDof,1);
nel=size(t,2);
ne=size(e,2);
values_1=zeros(5,1);
values_2=zeros(5,1);
for iel=1:nel
    nd=t(1:3,iel);
    a=p(:,nd(1));
    b=p(:,nd(2));
    c=p(:,nd(3));
    f=elF(a,b,c);
    F(nd)=F(nd)+f;
end
end       

function [A,F]=boundaryMatrix(A,F,p,e,boundary_method)
%% boundary condition
values_1=zeros(5,1);
values_2=zeros(5,1);
% loc_neuEdges_2=find(p(1,e(1,:))==1);   %原始
loc_neuEdges_2=find(p(1,e(1,:))==1&p(1,e(2,:))==1);   %修改
neuEdges_2=e(:,loc_neuEdges_2);
nneu_2=size(neuEdges_2,2);
for iedge=1:nneu_2                          %第二条边
    y1=p(2,neuEdges_2(1,iedge));
    y2=p(2,neuEdges_2(2,iedge));
    min_y=min(y1,y2);
    max_y=max(y1,y2);
    y=get_integration_points_1d(min_y,max_y);
    for k=1:5
          values_1(k)=-nboundaryFunction_2(1,y(k))*(y(k)-y2)/(y1-y2);
    end
    F(neuEdges_2(1,iedge))=F(neuEdges_2(1,iedge))+gaussint5(min_y,max_y,values_1);
    for k=1:5
          values_2(k)=-nboundaryFunction_2(1,y(k))*(y(k)-y1)/(y2-y1);
    end
    F(neuEdges_2(2,iedge))=F(neuEdges_2(2,iedge))+gaussint5(min_y,max_y,values_2);
end

% loc_neuEdges_4=find(p(1,e(1,:))==-1);                  %原始
loc_neuEdges_4=find(p(1,e(1,:))==-1&p(1,e(2,:))==-1);    %修改
neuEdges_4=e(:,loc_neuEdges_4);
nneu_4=size(neuEdges_4,2);
for iedge=1:nneu_4                     %第四个边
    y1=p(2,neuEdges_4(1,iedge));
    y2=p(2,neuEdges_4(2,iedge));
    min_y=min(y1,y2);
    max_y=max(y1,y2);
    y=get_integration_points_1d(min_y,max_y);
    for k=1:5
          values_1(k)=nboundaryFunction_2(-1,y(k))*(y(k)-y2)/(y1-y2);
    end
    F(neuEdges_4(1,iedge))=F(neuEdges_4(1,iedge))+gaussint5(min_y,max_y,values_1);
    for k=1:5
          values_2(k)=nboundaryFunction_2(-1,y(k))*(y(k)-y1)/(y2-y1);
    end
    F(neuEdges_4(2,iedge))=F(neuEdges_4(2,iedge))+gaussint5(min_y,max_y,values_2);  
end

loc_d=find(p(2,:)==-1|p(2,:)==1);
if boundary_method==1 
    x=p(1,loc_d);
    y=p(2,loc_d);
    A(loc_d,:)=0;
    A(loc_d,loc_d)=eye(size(loc_d,2));    %3*9变为3*3
    F(loc_d)=boundaryFunction(x,y)';
elseif boundary_method==2
    x=p(1,loc_d);
    y=p(2,loc_d);
    F=F-A(:,loc_d)*boundaryFunction(x,y)';
    A(loc_d,:)=[];
    A(:,loc_d)=[];
    F(loc_d)=[];
else
    x=p(1,loc_d);
    y=p(2,loc_d);
    A(loc_d,loc_d)=10^15*eye(size(loc_d,2));
    F(loc_d)=10^15*boundaryFunction(x,y)';
end
end

function u=solveAF(A,F,p,e,boundary_method)
%% solve A\F
if boundary_method==1||boundary_method==3
    u=A\F;
else
    loc_d=find(p(2,:)==-1|p(2,:)==1);
    uInt=A\F;
    sDof=size(p,2);
    u=zeros(sDof,1);
    locInt=1:sDof;
    locInt(loc_d)=[];
    u(locInt)=uInt;
    x=p(1,loc_d)';
    y=p(2,loc_d)';
    uBoundary=boundaryFunction(x,y);
    u(loc_d)=uBoundary;
end
end

function [L2_error,H1_error]=errorEstimate(u,p,t)
%% Error estimate
% L2=sqrt(/iint(u-u_h)^2dxdy)

L2_error=0;
H1_half=0;
for i=1:size(t,2)
    triangle=t(1:3,i);
    a=p(:,triangle(1));
    b=p(:,triangle(2));
    c=p(:,triangle(3));
    w=[a b c];
    [x,y]=get_integration_points(a,b,c);
    values_numerial=tri_basis(a,b,c,a,x,y).*u(triangle(1))+...
        tri_basis(a,b,c,b,x,y).*u(triangle(2))+...
        tri_basis(a,b,c,c,x,y).*u(triangle(3));
    values_true=exactSolution(x,y);
    values=(values_numerial-values_true).^2;
    L2_error=L2_error+triangle_quadrature(a,b,c,values);
    
    values_numerial_gradient=tri_basis_gradient(a,b,c,a,x,y).*u(triangle(1))+...
        tri_basis_gradient(a,b,c,b,x,y).*u(triangle(2))+...
        tri_basis_gradient(a,b,c,c,x,y).*u(triangle(3));
    exactValues_gradient=exactSolution_gradient(x,y);
    values_gradient=sum((values_numerial_gradient-exactValues_gradient).^2,2);
    H1_half=H1_half+triangle_quadrature(a,b,c,values_gradient);    
end
L2_error=sqrt(L2_error);
H1_error=sqrt(L2_error^2+H1_half);
end


function plotFigure(p,e,t,u)
%% Plot figure
x=p(1,:)';
y=p(2,:)';
z=exactSolution(x,y);
subplot(1,2,1);
pdemesh(p,e,t,u);
title('numerical solution');
subplot(1,2,2);
pdemesh(p,e,t,z);
title('exact solution');
end



function z=exactSolution(x,y)
%% Exact solution
% z=sin(pi.*x).*sin(pi.*y);
z=cos(pi.*x).*cos(pi.*y);
end

function z=exactSolution_gradient(x,y)
%% The gradient of exact solution
% z(:,1)=pi.*cos(pi.*x).*sin(pi.*y);
% z(:,2)=pi.*sin(pi.*x).*cos(pi.*y);
z(:,1)=-pi.*sin(pi.*x).*cos(pi.*y);
z(:,2)=-pi.*cos(pi.*x).*sin(pi.*y);
end

function z=boundaryFunction(x,y)
%% Boundary function
% z=sin(pi.*x).*sin(pi.*y);
z=cos(pi.*x).*cos(pi.*y);
end

function z=nboundaryFunction_1(x,y)
%% Boundary function
% z=sin(pi.*x).*sin(pi.*y);
z=pi*cos(pi.*x).*sin(pi.*y);
end

function z=nboundaryFunction_2(x,y)
%% Boundary function
% z=sin(pi.*x).*sin(pi.*y);
z=pi*sin(pi.*x).*cos(pi.*y);
end

function z=rightFunction(x,y)
% z=2*pi*pi*sin(pi.*x).*sin(pi.*y);
z=2*pi*pi*cos(pi.*x).*cos(pi.*y);
end

function k_f=elF(a,b,c)
%% Local right hands
w=[a b c];
[x,y]=get_integration_points(a,b,c);
for j=1:3
    basis_values_j=tri_basis(a,b,c,w(:,j),x,y);
    right_values=rightFunction(x,y);
    values=sum(basis_values_j.*right_values,2);
    k_f(j)=triangle_quadrature(a,b,c,values);
end
k_f=k_f';
end

function k_matrix=elK(a,b,c)
%% Generate local stiffness matrix
k_matrix=zeros(3,3);
w=[a b c];
[x,y]=get_integration_points(a,b,c);
for j=1:3
    gradient_values_j=tri_basis_gradient(a,b,c,w(:,j),x,y);
    for k=1:3
        gradient_values_k=tri_basis_gradient(a,b,c,w(:,k),x,y);
        gradient_values=sum(gradient_values_j.*gradient_values_k,2);
        k_matrix(j,k)=triangle_quadrature(a,b,c,gradient_values);
    end
end
end

function z=tri_basis(a,b,c,point,x,y)
%% Basis function
x1=a(1); y1=a(2); 
x2=b(1); y2=b(2);
x3=c(1); y3=c(2);
S=0.5*(x2*y3+x1*y2+x3*y1-x2*y1-x1*y3-x3*y2);
if point==a
    Sxy=0.5*(x2*y3+x*y2+x3*y-x2*y-x*y3-x3*y2);
elseif point==b
    Sxy=0.5*(x*y3+x1*y+x3*y1-x*y1-x1*y3-x3*y);
else
    Sxy=0.5*(x2*y+x1*y2+x*y1-x2*y1-x1*y-x*y2);
end
z=Sxy/S;
end

function z=tri_basis_gradient(a,b,c,point,x,y)
%% Gradient function
x1=a(1); y1=a(2); 
x2=b(1); y2=b(2);
x3=c(1); y3=c(2);
S=0.5*(x2*y3+x1*y2+x3*y1-x2*y1-x1*y3-x3*y2);
z=zeros(length(x),2);
w=[a b c];
num=0;
x_vector=zeros(2,1);
y_vector=zeros(2,1);
if point==a
    num=1;
elseif point==b
    num=2;
else
    num=3;
end
k=mod(num,3)+1;
x_vector(1)=w(1,k);
y_vector(1)=w(2,k);
k=k+1;
if(k>3)
    k=1;
end
x_vector(2)=w(1,k);
y_vector(2)=w(2,k);
z(:,1)=(y_vector(1)-y_vector(2))/S/2;
z(:,2)=(x_vector(2)-x_vector(1))/S/2;
end

function [x,y]=get_integration_points(a,b,c)
%% Generate intergration points
x_temp=[a(1);b(1);c(1)];
y_temp=[a(2);b(2);c(2)];
x=[x_temp;sum(x_temp)/3;sum(x_temp([1,2]))/2;sum(x_temp([2,3]))/2;sum(x_temp([1,3]))/2];
y=[y_temp;sum(y_temp)/3;sum(y_temp([1,2]))/2;sum(y_temp([2,3]))/2;sum(y_temp([1,3]))/2];
end

function x=get_integration_points_1d(a,b)
%% Generate intergration points(1D)
location=[0.046910077030668 0.230765344947158 0.5  0.769234655052842 0.953089922969332];
for i=1:5
x(i)=a+location(i)*(b-a);
end
end

function y=gaussint5(a,b,values)
weights=[0.236926885056189 0.478628670499366 0.568888888888889 0.478628670499366 0.236926885056189];
y=(weights*values)*((b-a)/2);
end

function y = triangle_quadrature(a,b,c,values)
%% Quadrature formula
weights = [1/20;1/20;1/20;9/20;2/15;2/15;2/15];
x1=a(1); y1=a(2); 
x2=b(1); y2=b(2);
x3=c(1); y3=c(2);
S=0.5*(x2*y3+x1*y2+x3*y1-x2*y1-x1*y3-x3*y2);
y=S*(weights'*values);
end

function [p,e,t]=uniform_trimesh_delaunay(x_start,x_end,y_start,y_end,n_x,n_y)
h_x=(x_end-x_start)/n_x;
h_y=(y_end-y_start)/n_y;

[X,Y] = meshgrid(x_start:h_x:x_end,y_start:h_y:y_end);
length=(n_x+1)*(n_y+1);
x=reshape(X,length,1);
y=reshape(Y,length,1);
p=[x,y]';
t=delaunay(x,y)';
t=[t;ones(1,size(t,2))];
loc_1=find(p(1,:)==x_start);
loc_2=find(p(1,:)==x_end);
loc_3=find(p(2,:)==y_start);
loc_4=find(p(2,:)==y_end);
e=[loc_1(1:n_y),loc_2(1:n_y),loc_3(1:n_x),loc_4(1:n_x);loc_1(2:n_y+1),loc_2(2:n_y+1),loc_3(2:n_x+1),loc_4(2:n_x+1)];
end

function [L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(pde,N)
%% Compute the convergence order
L2_error_array=zeros(N,1);
H1_error_array=zeros(N,1);
num_freedom_array=zeros(N,1);
L2_convergenceOrder=zeros(N-1,1);
H1_convergenceOrder=zeros(N-1,1);
N_vector=2.^(2:N+1);

for i=1:N
    [num_freedom,L2_error,H1_error]=run_main(N_vector(i),pde);
    num_freedom_array(i)=num_freedom;
    L2_error_array(i)=L2_error;
    H1_error_array(i)=H1_error;
end
for i=1:N-1
    L2_convergenceOrder(i)=log2(L2_error_array(i)/L2_error_array(i+1));
    H1_convergenceOrder(i)=log2(H1_error_array(i)/H1_error_array(i+1));
end
figure(2)
loglog(num_freedom_array,H1_error_array,'-o',num_freedom_array,exp(-1*log(sqrt(num_freedom_array))+1.9),'',...
    num_freedom_array,L2_error_array,'-*',num_freedom_array,exp(-2*log(sqrt(num_freedom_array))+1.5),'-.')
xlabel('degree of freedom in log scale');
ylabel('error in log scale');
hleg1=legend('H1-norm','slope of order 1','L2-norm','slope of order 2','Location','SouthWest');
title('Convergence Order')
%axis([1e2,1e5,1e-4,1e1])
end



