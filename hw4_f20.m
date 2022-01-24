function hw4_f20
clear; clc; close all; format compact; format short e
%
Problem=input('Enter desired problem (1-3): ');
%
% ----------------------------------------------------------------
%
if Problem==1
    disp('Problem 1 - 4.7')
    E=30e6;
    I=200;
    L=240; % length in inches
    K=k_global(E,I,L)
    % Nodal displacements and rotations:
    F=[-500;0;0]
    Kp=K([1,2,4],[1,2,4])
    d=Kp\F;
    v1=d(1)
    phi1=d(2)
    phi2=d(3)
    % Reactions
    Reactions=K*[v1;phi1;0;phi2;0;0]
    % Element loads
    k=k_element(E,I,L);
    f1=k*[v1;phi1;0;phi2]
    f2=k*[0;phi2;0;0]
    %
    % ----------------------------------------------------------------
    %
elseif Problem==2
    disp('Problem 2')
    E=29e6;
    I=200;
    L=15*12; % length in inches
    K=k_global(E,I,L)
    % Nodal displacements and rotations:
    w=1000/12;  % lb/in
    F=[-w*L;0;w*L^2/12]
    Kp=K([3,4,6],[3,4,6])
    d=Kp\F;
    v2=d(1)
    phi2=d(2)
    phi3=d(3)
    % Reactions
    Reactions=K*[0;0;v2;phi2;0;phi3]-[-w*L/2;-w*L^2/12;-w*L;0;-w*L/2;w*L^2/12]
    % Element loads
    k=k_element(E,I,L);
    f1=k*[0;0;v2;phi2]-[-w*L/2;-w*L^2/12;-w*L/2;w*L^2/12]
    f2=k*[v2;phi2;0;phi3]-[-w*L/2;-w*L^2/12;-w*L/2;w*L^2/12]
    %
    % ----------------------------------------------------------------
    %
elseif Problem==3
    disp('Problem 3')
    E=29e6;
    I=150;
    L=120; % length in inches
    K=k_global(E,I,L)
    % Nodal displacements and rotations:
    w=2000/12; % lb/in
    % work equivalent nodal loads
    F1y=-3*w*L/20;
    M1=-w*L^2/30;
    F2y=-w*L;
    M2=-w*L^2/15;
    F3y=-17*w*L/20;
    M3=2*w*L^2/15;
    F=[M2;F3y;M3]
    Kp=K([4:6],[4:6])
    d=Kp\F;
    phi2=d(1)
    v3=d(2)
    phi3=d(3)
    % Reactions
    Reactions=K*[0;0;0;phi2;v3;phi3]-[F1y; M1; F2y; M2; F3y; M3]
    % Element loads
    k=k_element(E,I,L);
    f1=k*[0;0;0;phi2]-[-3*w*L/20;-w*L^2/30;-7*w*L/20;w*L^2/20]
    f2=k*[0;phi2;v3;phi3]-[-13*w*L/20;-7*w*L^2/60;-17*w*L/20;2*w*L^2/15]
    %
    % ----------------------------------------------------------------
    %
else
    disp('Problem 4, part b')
    %
    E=29e6;
    I=200;
    L=15*12;  % length in inches
    % FEA solution
    % element 1
    x=linspace(0,L,101);
    d=[0;0;-1.2569;-3.4914e-03];
    for i=1:101
        xp(i)=x(i);
        N=interp(x(i),L);
        v(i)=N*d;
    end
    plot(xp,v,'b')
    hold on
    % element 2
    d=[-1.2569;-3.4914e-03; 0; .013966];
    for i=1:101
        xp(i)=x(i)+L;
        N=interp(x(i),L);
        v(i)=N*d;
    end
    plot(xp,v,'g')
    % exact
    l=30*12; % length in inches
    w0=1000/12;
    xe=linspace(0,l,101);
    ve=-((w0*l^4)/(E*I))*((1/16)*(xe/l).^2-(5/48)*(xe/l).^3+(1/24)*(xe/l).^4);
    plot(xe,ve,'r')
    legend('element 1','element 2','exact')
    xlabel('x')
    ylabel('v(x)')
    title('Problem 4')
end
%
%-----------------------------------------------------------------
%
function K=k_global(E,I,L)
%
k1=zeros(6,6);
k2=zeros(6,6);
k1(1:4,1:4)=k_element(E,I,L);
k2(3:6,3:6)=k_element(E,I,L);
K=k1+k2;
%
%-----------------------------------------------------------------
%
function k=k_element(E,I,L)
%
k=(E*I/L^3)*[12,6*L,-12,6*L;
    6*L,4*L^2,-6*L,2*L^2;
    -12,-6*L,12,-6*L;
    6*L,2*L^2,-6*L,4*L^2]
%
%-----------------------------------------------------------------
%
function N=interp(x,L)
%
N1=(2*x.^3-3*x.^2*L+L^3)/L^3;
N2=(x.^3*L-2*x.^2*L^2+x*L^3)/L^3;
N3=(-2*x.^3+3*x.^2*L)/L^3;
N4=(x.^3*L-x.^2*L^2)/L^3;
N=[N1, N2, N3, N4];
