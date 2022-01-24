function example_6_1
%
clc; clear all; close all
%
E=30e6;
nu=0.25;
t=1;
plane_stress=true; % Note: use plane_stress=false for plane strain problems
x=[0 2 0];
y=[-1 0 1];
% stiffness matrix
k=k_cst(E,nu,t,x,y,plane_stress)
% stresses
u=[0 .0012 0];
v=[.0025 0 .0025];
sigma=sigma_cst(E,nu,t,x,y,plane_stress,u,v)
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
