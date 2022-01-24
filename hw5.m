% hw 5 1) 6.3c, 2) 6.11b, 3) 6.13
clc; clear; format short
%
% 1
%
disp('Problem 1')
E=30e6;
nu=0.25;
t=1;
xi=0; yi=0;
xj=2; yj=0;
xm=0; ym=1;
betai=yj-ym;  betaj=ym-yi;  betam=yi-yj;
gammai=xm-xj; gammaj=xi-xm; gammam=xj-xi;
A=.5*det([1 xi yi; 1 xj yj; 1 xm ym]);
B=(1/(2*A))*[betai 0 betaj 0 betam 0;
    0 gammai 0 gammaj 0 gammam;
    gammai betai gammaj betaj gammam betam]
D=(E/(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]
k=t*A*B'*D*B
%
% 2
%
%   SOLVED BY HAND - SEE NOTES
%
% 3
%
disp('Problem 3')
clear; format short e
E=30e6;
nu=0.3;
t=1;
D=(E/(1-nu^2))*[1 nu 0; nu 1 0; 0 0 (1-nu)/2]
% From section 6.5, page 270, with u1=v1=u2=u4, we use rows 5-8 and 
%  columns 5-8 of the K matrix, which yields:
K=(375000/0.91)*[48  0  -28  14;
                  0  87  12 -80;
                -28  12  48 -26;
                 14 -80 -26  87];
% the net shear load, 10,000 lb is split equally between nodes 3 & 4:
F=[0; -5000; 0; -5000];
% solver for nodal displacements
d=K\F
% B matrices
% element 1, from p. 366
B1=(1/200)*[0 0 10 0 -10 0; 0 -20 0 0 0 20; -20 0 0 10 20 -10];
% element 2, from p. 366
B2=(1/200)*[-10 0 10 0 0 0; 0 0 0 -20 0 20; 0 -10 -20 10 20 0];
% stresses - element 1
disp('Element 1 stresses')
sigma=D*B1*[0; 0; d(1); d(2); 0; 0]
sx=sigma(1); sy=sigma(2); txy=sigma(3);
sigma_1=(sx+sy)/2+sqrt((((sx-sy)/2))^2+txy^2)
sigma_2=(sx+sy)/2-sqrt((((sx-sy)/2))^2+txy^2)
% stresses - element 2
disp('Element 2 stresses')
sigma=D*B2*[0; 0; d(3); d(4); d(1); d(2)]
sx=sigma(1); sy=sigma(2); txy=sigma(3);
sigma_1=(sx+sy)/2+sqrt((((sx-sy)/2))^2+txy^2)
sigma_2=(sx+sy)/2-sqrt((((sx-sy)/2))^2+txy^2)
                        
