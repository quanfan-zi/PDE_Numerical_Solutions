clear;clc;
format long
% 主参数
k=10; %k表示循环次数，n表示基函数个数 
Con1=zeros(1,k);
Con2=zeros(1,k);  %Coni储存矩阵条件数
Error1=zeros(1,k);
Error2=zeros(1,k); %Error储存L2误差

for n=1:k
% 系数矩阵
A1=zeros(n,n);
A2=zeros(n,n);
for i=1:n
    for j=1:n
        if i==j
            A1(i,j)=1/2*((i*pi)^2-1);
        end
        A2(i,j)=i*j/(i+j-1)-2*i*j/(i+j)-1+(i*j+i+j)/(i+j+1)+2/(i+j+2)-1/(i+j+3);
    end
end
Con1(n)=cond(A1,2);
Con2(n)=cond(A2,2);

% 右端项
b1=zeros(n,1);
b2=zeros(n,1);
for i=1:n
    b1(i)=-cos(i*pi)/(i*pi);
    b2(i)=1/(i+2)-1/(i+3);
end
w1=A1\b1;
w2=A2\b2;

% 计算误差
x=0:0.01:1;   %取101个节点
r1=zeros(1,101);
r2=zeros(1,101);
for i=1:n
    r1=r1+w1(i)*Base1(x,i);
    r2=r2+w2(i)*Base2(x,i);
end
r1=RealRes(x)-r1;
r2=RealRes(x)-r2;

Error1(n)=(ComSimpson(0,1,r1.^2))^0.5;
Error2(n)=(ComSimpson(0,1,r2.^2))^0.5;
end

%画图
N=1:k;

subplot(2,2,1);
plot(Error1,'b','LineWidth',1);
title('\phi(x)与真解的L^2误差')
subplot(2,2,3);
semilogy(Error1,'b','Linewidth',1)
title('\phi(x)与真解的L^2误差(对数图)')
subplot(2,2,2);
plot(Error2,'r','Linewidth',1);
title('\psi(x)与真解的L^2误差')
subplot(2,2,4);
semilogy(Error2,'r','Linewidth',1)
title('\psi(x)与真解的L^2误差(对数图)')

function I = ComSimpson(a, b, f)
% Simpson公式

N=length(f)-1;
h=(b-a)/N;

y=2*ones(1, N+1);
y(2:2:N)=4;
y(1)=1;
y(N+1)=1;

I=f*y'*h/3;
end

function y=RealRes(x)
% 真解
y=sin(x)/sin(1)-x;
end

function y=Base1(x,i)
% 第一类基函数
y=sin(i*pi*x);
end

function y=Base2(x,i)
% 第二类基函数
y=(1-x).*((x).^i);
end









