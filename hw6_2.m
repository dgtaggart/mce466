function example_6_2
%
clc; clear all; close all
%
E=30e6;
nu=0.3;
t=1;
plane_stress=true; % Note: use plane_stress=false for plane strain problems
for element=1:2
    if element==1
        x=[0 20  0];
        y=[0 10 10];
        k1=k_cst(E,nu,t,x,y,plane_stress);
        % rearrange rows and columns
        k1p=[k1(:,1:2),k1(:,5:6),k1(:,3:4)];
        k1p=[k1p(1:2,:);k1p(5:6,:);k1p(3:4,:)]
    else
        x=[0 20 20];
        y=[0  0 10];
        k2=k_cst(E,nu,t,x,y,plane_stress);
        % rearrange rows and columns
        k2p=[k2(:,1:2),k2(:,5:6),k2(:,3:4)];
        k2p=[k2p(1:2,:);k2p(5:6,:);k2p(3:4,:)]
    end
end
K=zeros(8,8);
K(1:6,1:6)=K(1:6,1:6)+k1p;
K([1,2,5:8],[1,2,5:8])=K([1,2,5:8],[1,2,5:8])+k2p;
P_known=[5000; 0; 5000; 0];
d_unknown=K(5:8,5:8)\P_known
d=[0;0;0;0;d_unknown]
for element=1:2
    if element==1
        x=[0 20  0];
        y=[0 10 10];
        u=[d(1) d(5) d(3)];
        v=[d(2) d(6) d(4)];
        sigma=sigma_cst(E,nu,t,x,y,plane_stress,u,v)
    else
        x=[0 20 20];
        y=[0  0 10];
        u=[d(1) d(7) d(5)];
        v=[d(2) d(8) d(6)];
        sigma=sigma_cst(E,nu,t,x,y,plane_stress,u,v)
    end
end

%
function k=k_cst(E,nu,t,x,y,plane_stress)
%
if plane_stress==true
    D=E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    D=E/((1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2];
end
%
A=.5*det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
beta=[y(2)-y(3) y(3)-y(1) y(1)-y(2)];
gamma=[x(3)-x(2) x(1)-x(3) x(2)-x(1)];
B=[ beta(1)    0      beta(2)    0       beta(3)    0;
    0      gamma(1)   0      gamma(2)    0      gamma(3);
    gamma(1) beta(1)  gamma(2)  beta(2)  gamma(3) beta(3)]/(2*A);
k=t*A*B'*D*B;
%
function sigma=sigma_cst(E,nu,t,x,y,plane_stress,u,v)
%
if plane_stress==true
    D=E/(1-nu^2)*[1 nu 0; nu 1 0; 0 0 (1-nu)/2];
else
    D=E/((1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0; 0 0 (1-2*nu)/2];
end
%
A=.5*det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
beta=[y(2)-y(3) y(3)-y(1) y(1)-y(2)];
gamma=[x(3)-x(2) x(1)-x(3) x(2)-x(1)];
B=[ beta(1)    0      beta(2)    0       beta(3)    0;
    0      gamma(1)   0      gamma(2)    0      gamma(3);
    gamma(1) beta(1)  gamma(2)  beta(2)  gamma(3) beta(3)]/(2*A);

d=[u(1); v(1); u(2); v(2); u(3); v(3)];
sigma=D*B*d;
