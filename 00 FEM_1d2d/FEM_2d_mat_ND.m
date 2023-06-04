function FEM_2d_mat_ND
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


N=9;
[L2_convergenceOrder,H1_convergenceOrder]=convergenceOrder(pde,N)

 %M=20;
% 
 %[sDof,L2_error,H1_error]=run_main(M,pde)

toc
end

function [sDof,L2_error,H1_error]=run_main(M,pde)
%% 

boundary_method=1;  % 1  2  3
domainLenght=pde.end_point_x-pde.start_point_x;
hmax=domainLenght/M;

% Generate mesh
[p,e,t]=meshGeneration(pde,M,hmax);
sDof=size(p,2);

% Generate stiffness matrix
A=stiffnessMatrix(p,t);

% Generate right hands
F=rightHands(p,t);

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
% g=[ 2     2     2     2
%     -1     1     1    -1
%      1     1    -1    -1
%     -1    -1     1     1
%     -1     1     1    -1
%      1     1     1     1
%      0     0     0     0];
% [p,e,t]=initmesh(g,'Hmax',hmax);
[p,e,t]=uniform_trimesh_delaunay(pde.start_point_x,pde.end_point_x,pde.start_point_y,pde.end_point_y,M,M);
end

function A=stiffnessMatrix(p,t)
%% Generate stiffness matrix
sDof=size(p,2);
nel=size(t,2);
nd=t(1:3,:);                      %3*nel (nel是三角形个数)
a=p(:,nd(1,:));                   %2*nel
b=p(:,nd(2,:));
c=p(:,nd(3,:));
x1=a(1,:); y1=a(2,:);
x2=b(1,:); y2=b(2,:);
x3=c(1,:); y3=c(2,:);
S=0.5.*(x2.*y3+x1.*y2+x3.*y1-x2.*y1-x1.*y3-x3.*y2);       %1*nel
S_sqrt=2.*sqrt(S);
Dx=[(y2-y3)./S_sqrt;(y3-y1)./S_sqrt;(y1-y2)./S_sqrt];     %3*nel 三个基底关于x求导
Dy=[(x3-x2)./S_sqrt;(x1-x3)./S_sqrt;(x2-x1)./S_sqrt];

DiDj_x=[Dx([1 2 3],:).*Dx([1 2 3],:);Dx([1 2 3],:).*Dx([2 3 1],:);Dx([1 2 3],:).*Dx([3 1 2],:)];   %9*nel
DiDj_y=[Dy([1 2 3],:).*Dy([1 2 3],:);Dy([1 2 3],:).*Dy([2 3 1],:);Dy([1 2 3],:).*Dy([3 1 2],:)];
DiDj=reshape((DiDj_x+DiDj_y),1,9*nel);
ii=reshape([nd;nd;nd],1,9*nel);
jj=reshape([nd([1 2 3],:);nd([2 3 1],:);nd([3 1 2],:)],1,9*nel);
A=sparse(ii,jj,DiDj,sDof,sDof);

% clear DiDj DiDj_x DiDj_y ii jj;                       %内存很大，但子函数

end

function F=rightHands(p,t)
%% Generate right hands
sDof=size(p,2);
nel=size(t,2);
nd=t(1:3,:);                      %3*nel
a=p(:,nd(1,:));
b=p(:,nd(2,:));
c=p(:,nd(3,:));
[x,y]=get_integration_points(a,b,c);             %x,y  7*nel
right_values=rightFunction(x,y);
basis_values_a=tri_basis(a,b,c,a,x,y);
basis_values_b=tri_basis(a,b,c,b,x,y);
basis_values_c=tri_basis(a,b,c,c,x,y);

values_a=triangle_quadrature(a,b,c,right_values.*basis_values_a)';           %1*nel
values_b=triangle_quadrature(a,b,c,right_values.*basis_values_b)';           %1*nel
values_c=triangle_quadrature(a,b,c,right_values.*basis_values_c)';
values=reshape([values_a;values_b;values_c],1,3*nel);                        %3*nel
ii=reshape(nd,1,3*nel);
F=accumarray(ii',values',[sDof,1]);
clear ii values values_a values_b values_c basis_values_a basis_values_b basis_values_c right_values;
end

function [A,F]=boundaryMatrix(A,F,p,e,boundary_method)
%% boundary condition
% loc_n_2=find(p(1,e(1,:))==1);    %原始
% loc_n_4=find(p(1,e(1,:))==-1);   %原始
loc_n_2=find(p(1,e(1,:))==1&p(1,e(2,:))==1);    %修改
loc_n_4=find(p(1,e(1,:))==-1&p(1,e(2,:))==-1);  %修改
loc_d=find(p(2,:)==-1|p(2,:)==1);
%ed=e(:,loc_d);
en_2=e(:,loc_n_2);
en_4=e(:,loc_n_4);
%%%%%%%%%%%%%%%%%% Neumman boundary condition %%%%%%%%%%%%%%%%%%%%%%%

    y1=p(2,en_2(1,:))';
    y2=p(2,en_2(2,:))';
    min_y=min(y1,y2);
    max_y=max(y1,y2);
    y=get_integration_points_1d(min_y,max_y);
    values_1=-nboundaryFunction_2(ones(size(y)),y).*(y-repmat(y2,1,size(y,2)))./(repmat(y1,1,size(y,2))-repmat(y2,1,size(y,2)));
    values_2=-nboundaryFunction_2(ones(size(y)),y).*(y-repmat(y1,1,size(y,2)))./(repmat(y2,1,size(y,2))-repmat(y1,1,size(y,2)));
    F(en_2(1,:))=F(en_2(1,:))+gaussint5(min_y,max_y,values_1);
    F(en_2(2,:))=F(en_2(2,:))+gaussint5(min_y,max_y,values_2);
    
    y1=p(2,en_4(1,:))';
    y2=p(2,en_4(2,:))';
    min_y=min(y1,y2);
    max_y=max(y1,y2);
    y=get_integration_points_1d(min_y,max_y);
    values_1=nboundaryFunction_2(-ones(size(y)),y).*(y-repmat(y2,1,size(y,2)))./(repmat(y1,1,size(y,2))-repmat(y2,1,size(y,2)));
    values_2=nboundaryFunction_2(-ones(size(y)),y).*(y-repmat(y1,1,size(y,2)))./(repmat(y2,1,size(y,2))-repmat(y1,1,size(y,2)));
    F(en_4(1,:))=F(en_4(1,:))+gaussint5(min_y,max_y,values_1);
    F(en_4(2,:))=F(en_4(2,:))+gaussint5(min_y,max_y,values_2);
%%%%%%%%%%%%%%%%%%%%%%%%%% Direchlet boundary condition %%%%%%%%%%%%%%%%%
if boundary_method==1 
    x=p(1,loc_d);
    y=p(2,loc_d);
    A(loc_d,:)=0;
    A(loc_d,loc_d)=eye(size(loc_d,2));
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
triangle=t(1:3,:);
a=p(:,triangle(1,:));
b=p(:,triangle(2,:));
c=p(:,triangle(3,:));
x1=a(1,:); y1=a(2,:);
x2=b(1,:); y2=b(2,:);
x3=c(1,:); y3=c(2,:);
S=0.5.*(x2.*y3+x1.*y2+x3.*y1-x2.*y1-x1.*y3-x3.*y2);
Dx=[(y2-y3)./(2.*S);(y3-y1)./(2.*S);(y1-y2)./(2.*S)];
Dy=[(x3-x2)./(2.*S);(x1-x3)./(2.*S);(x2-x1)./(2.*S)];
[x,y]=get_integration_points(a,b,c);
col_x=size(x,2);
values_numerial=tri_basis(a,b,c,a,x,y).*repmat(u(triangle(1,:)),1,col_x)+...
        tri_basis(a,b,c,b,x,y).*repmat(u(triangle(2,:)),1,col_x)+...
        tri_basis(a,b,c,c,x,y).*repmat(u(triangle(3,:)),1,col_x);
values_true=exactSolution(x,y);
values=(values_numerial-values_true).^2;
L2_error=sum(triangle_quadrature(a,b,c,values));

values_numerial_x=repmat(sum(Dx.*u(triangle),1)',1,col_x);
values_numerial_y=repmat(sum(Dy.*u(triangle),1)',1,col_x);
[values_exat_x,values_exat_y]=exactSolution_gradient(x,y);
values_gradient=(values_numerial_x-values_exat_x).^2+(values_numerial_y-values_exat_y).^2;
H1_half=sum(triangle_quadrature(a,b,c,values_gradient));

L2_error=sqrt(L2_error);
H1_error=sqrt(L2_error^2+H1_half);

end


function plotFigure(p,e,t,u)
%% Plot figure
figure(1)
pdemesh(p,e,t)
figure(2)
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

function [z_x,z_y]=exactSolution_gradient(x,y)
%% The gradient of exact solution
% z(:,1)=pi.*cos(pi.*x).*sin(pi.*y);
% z(:,2)=pi.*sin(pi.*x).*cos(pi.*y);
z_x=-pi.*sin(pi.*x).*cos(pi.*y);
z_y=-pi.*cos(pi.*x).*sin(pi.*y);
end

function z=boundaryFunction(x,y)
%% Boundary function
% z=sin(pi.*x).*sin(pi.*y);
z=cos(pi.*x).*cos(pi.*y);
end

function z=rightFunction(x,y)
% z=2*pi*pi*sin(pi.*x).*sin(pi.*y);
z=2*pi*pi*cos(pi.*x).*cos(pi.*y);
end


function z=tri_basis(a,b,c,point,x,y)
%% Basis function
col_x=size(x,2);
x1=repmat(a(1,:)',1,col_x); y1=repmat(a(2,:)',1,col_x); 
x2=repmat(b(1,:)',1,col_x); y2=repmat(b(2,:)',1,col_x);
x3=repmat(c(1,:)',1,col_x); y3=repmat(c(2,:)',1,col_x);
S=0.5.*(x2.*y3+x1.*y2+x3.*y1-x2.*y1-x1.*y3-x3.*y2);

if point==a
    Sxy=0.5.*(x2.*y3+x.*y2+x3.*y-x2.*y-x.*y3-x3.*y2);
elseif point==b
    Sxy=0.5.*(x.*y3+x1.*y+x3.*y1-x.*y1-x1.*y3-x3.*y);
else
    Sxy=0.5.*(x2.*y+x1.*y2+x.*y1-x2.*y1-x1.*y-x.*y2);
end
z=Sxy./S;
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

function x=get_integration_points_1d(a,b)
%% Generate intergration points(1D)
location=[0.046910077030668 0.230765344947158 0.5  0.769234655052842 0.953089922969332];
h=b-a;
x=repmat(a,1,5)+h*location;
end

function y=gaussint5(a,b,values)
weights=[0.236926885056189;0.478628670499366;0.568888888888889;0.478628670499366;0.236926885056189];
y=(values*weights).*((b-a)/2);
end


function [x,y]=get_integration_points(a,b,c)
%% Generate intergration points
x_temp=[a(1,:);b(1,:);c(1,:)];
y_temp=[a(2,:);b(2,:);c(2,:)];
x=[x_temp;sum(x_temp,1)/3;sum(x_temp([1,2],:),1)/2;sum(x_temp([2,3],:),1)/2;sum(x_temp([1,3],:),1)/2]';
y=[y_temp;sum(y_temp,1)/3;sum(y_temp([1,2],:),1)/2;sum(y_temp([2,3],:),1)/2;sum(y_temp([1,3],:),1)/2]';
end

function y = triangle_quadrature(a,b,c,values)
%% Quadrature formula
weights = [1/20;1/20;1/20;9/20;2/15;2/15;2/15];
x1=a(1,:); y1=a(2,:); 
x2=b(1,:); y2=b(2,:);
x3=c(1,:); y3=c(2,:);
S=0.5.*(x2.*y3+x1.*y2+x3.*y1-x2.*y1-x1.*y3-x3.*y2);
y=(values*weights).*S';
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
