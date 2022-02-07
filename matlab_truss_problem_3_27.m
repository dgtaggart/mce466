function truss_2D
%
%----------------------------------------------------------------------------%
% Static 2-d truss solution
%
%   adapted from the text:
%     The Finite Element Method Using MATLAB
%           by Young W. Kwon and Hyochoong Bang
%
%   adapted by D. G. Taggart for MCE 466
%
% Variable descriptions
%   k = element stiffness matrix
%   kk = system stiffness matrix
%   ff = system force vector
%   index = a vector containing system dofs associated with each element
%   gcoord = global coordinate matrix
%   disp = nodal displacement vector
%   elforce = element force vector
%   eldisp = element nodal displacement
%   stress = stress vector for every element
%   elprop = element property matrix
%   nodes = nodal connectivity matrix for each element
%   bcdof = a vector containing dofs associated with boundary conditions
%   bcval = a vector containing boundary condition values associated with
%           the dofs in 'bcdof'
%----------------------------------------------------------------------------%
clear; clc; format long; close all
%-------------------- Begin model definition --------------------
%
%---------------------------
% joint definitions (x,y,BC,Fx,Fy)
%    where BC=0 (u&v-free)
%    where BC=1 (u=v=0)
%    where BC=2 (u=0,v-free)
%    where BC=3 (u-free,v=0)
%---------------------------
node_def=[ -12*30,     0,   1,      0,        0;
            12*30,     0,   1,      0,        0;
                0, 12*40,   0,   5000,   -20000;
                0, 12*60,   1,      0,        0];
%---------------------------
% member definitions (Node 1, Node 2, Area)
%---------------------------
elements=[1, 3, 3;
          2, 3, 3;
          3, 4, 3];
%---------------------------
% material property (Young's modulus)
%---------------------------
E=30e6;  
%
%-------------------- End model definition --------------------
%
%---------------------------
%  plot truss geometry
%---------------------------
plot_truss(node_def,elements)
%---------------------------
%  control input data
%---------------------------
nsiz=size(node_def);
esiz=size(elements);
nel=esiz(1);             % number of element
nnode=nsiz(1);           % number of nodes
nnel=2;                  % number of nodes per element
ndof=2;                  % number of dofs per node
sdof=nnode*ndof;         % total system dofs
%---------------------------
%  nodal coordinates
%---------------------------
gcoord=node_def(:,1:2);   % (x,y) coordinates of nodes
%------------------------------------------
%  material and geometric properties
%------------------------------------------
for iel=1:nel
    elprop(iel,1)=E;                % Young's modulus
    elprop(iel,2)=elements(iel,3);  % member cross-section
end
%-----------------------------
%  nodal connectivity
%-----------------------------
nodes=elements(:,1:2);
%-----------------------------
%  applied constraints
%-----------------------------
count=1;
for i=1:nnode
    if node_def(i,3)==1
        bcdof(1,count)=2*i-1;
        bcdof(1,count+1)=2*i;
        count=count+2;
    elseif node_def(i,3)==2
        bcdof(1,count)=2*i-1;
        count=count+1;
    elseif node_def(i,3)==3
        bcdof(1,count)=2*i;
        count=count+1;
    else
        continue
    end
end
bcval=zeros(1,count-1);
%----------------------------
%  initialization to zero
%----------------------------
ff=zeros(sdof,1);              % system force vector
kk=zeros(sdof,sdof);           % system stiffness matrix
index=zeros(nnel*ndof,1);      % index vector
elforce=zeros(nnel*ndof,1);    % element force vector
eldisp=zeros(nnel*ndof,1);     % element nodal displacement vector
k=zeros(nnel*ndof,nnel*ndof);  % element stiffness matrix
stress=zeros(nel,1);           % stress vector for every element
%-----------------------------
%  applied nodal force
%-----------------------------
for i=1:nnode
    idof=2*i-1;
    ff(idof)=node_def(i,4);
    ff(idof+1)=node_def(i,5);
end
%--------------------------
%  loop for elements
%--------------------------
for iel=1:nel             % loop for the total number of elements
    nd(1)=nodes(iel,1);   % 1st connected node for the (iel)-th element
    nd(2)=nodes(iel,2);   % 2nd connected node for the (iel)-th element
    x1=gcoord(nd(1),1);   % x coordinate of 1st node
    y1=gcoord(nd(1),2);   % y coordinate of 1st node
    x2=gcoord(nd(2),1);   % x coordinate of 1st node
    y2=gcoord(nd(2),2);   % y coordinate of 2nd node
    leng=sqrt((x2-x1)^2+(y2-y1)^2);   % element length
    if (x2-x1)==0
        beta=2*atan(1);   % angle between local and global axes
    else
        beta=atan((y2-y1)/(x2-x1));
    end
    el=elprop(iel,1);                  % extract elastic modulus
    area=elprop(iel,2);                % extract cross-sectional area
    index=feeldof(nd,nnel,ndof);       % extract system dofs for the element
    k=fetruss2(el,leng,area,0,beta,1); % compute element matrix
    kk=feasmbl1(kk,k,index);           % assemble into system matrix
end
%---------------------------------------------------
%  apply constraints and solve the matrix
%---------------------------------------------------
[kk,ff]=feaplyc2(kk,ff,bcdof,bcval);  % apply the boundary conditions
disp=kk\ff;   % solve the matrix equation to find nodal displacements
%--------------------------------------------------
%  post computation for stress calculation
%--------------------------------------------------
for iel=1:nel         % loop for the total number of elements
    nd(1)=nodes(iel,1);   % 1st connected node for the (iel)-th element
    nd(2)=nodes(iel,2);   % 2nd connected node for the (iel)-th element
    x1=gcoord(nd(1),1); y1=gcoord(nd(1),2);  % coordinate of 1st node
    x2=gcoord(nd(2),1); y2=gcoord(nd(2),2);  % coordinate of 2nd node
    leng=sqrt((x2-x1)^2+(y2-y1)^2);  % element length
    if (x2-x1)==0
        beta=2*atan(1);       % angle between local and global axes
    else
        beta=atan((y2-y1)/(x2-x1));
    end
    el=elprop(iel,1);               % extract elastic modulus
    area=elprop(iel,2);             % extract cross-sectional area
    index=feeldof(nd,nnel,ndof);    % extract system dofs for the element
    k=fetruss2(el,leng,area,0,beta,1); % compute element matrix
    for i=1:(nnel*ndof)             % extract displacements associated with
        eldisp(i)=disp(index(i));   % (iel)-th element
    end
    c=cos(beta); s=sin(beta);
    stress(iel)=(E/leng)*[-c -s c s]*eldisp;    
end
%----------------------------
% print fem solutions
%----------------------------
num=1:1:sdof;
displ=[num' disp]          % print displacements
numm=1:1:nel;
stresses=[numm' stress]    % print stresses
%
%--------------------------------------------------------------------
%
function [k,m]=fetruss2(el,leng,area,rho,beta,ipt)
%     Stiffness and mass matrices for the 2-d truss element
%     nodal dof {u_1 v_1 u_2 v_2}
%
%  Synopsis:
%     [k,m]=fetruss2(el,leng,area,rho,beta,ipt)
%
%  Variable Description:
%     k - element stiffness matrix (size of 4x4)
%     m - element mass matrix (size of 4x4)
%     el - elastic modulus
%     leng - element length
%     area - area of truss cross-section
%     rho - mass density (mass per unit volume)
%     beta - angle between the local and global axes         
%            positive if the local axis is in the ccw direction from
%            the global axis
%     ipt = 1 - consistent mass matrix
%         = 2 - lumped mass matrix
%--------------------------------------------------------------------------
% stiffness matrix
c=cos(beta); s=sin(beta);
k= (area*el/leng)*[ c*c   c*s  -c*c  -c*s;...
    c*s   s*s  -c*s  -s*s;...
    -c*c  -c*s   c*c   c*s;...
    -c*s  -s*s   c*s   s*s];
% consistent mass matrix
if ipt==1
    m=(rho*area*leng/6)*[ 2*c*c+2*s*s  0  c*c+s*s  0;...
        0  2*c*c+2*s*s  0  c*c+s*s;...
        c*c+s*s  0  2*c*c+2*s*s  0;...
        0  c*c+s*s  0  2*c*c+2*s*s];
    % lumped mass matrix
else
    m=(rho*area*leng/2)*[ c*c+s*s  0  0  0;...
        0  c*c+s*s  0  0;...
        0  0  c*c+s*s  0;...
        0  0  0  c*c+s*s];
end
%
%--------------------------------------------------------------------
%
function [kk]=feasmbl1(kk,k,index)
%  Purpose:
%     Assembly of element matrices into the system matrix
%
%  Synopsis:
%     [kk]=feasmbl1(kk,k,index)
%
%  Variable Description:
%     kk - system matrix
%     k  - element matri
%     index - d.o.f. vector associated with an element
%-----------------------------------------------------------
edof = length(index);
for i=1:edof
    ii=index(i);
    for j=1:edof
        jj=index(j);
        kk(ii,jj)=kk(ii,jj)+k(i,j);
    end
end
%
%--------------------------------------------------------------------
%
function [index]=feeldof(nd,nnel,ndof)
%  Purpose:
%     Compute system dofs associated with each element
%
%  Synopsis:
%     [index]=feeldof(nd,nnel,ndof)
%
%  Variable Description:
%     index - system dof vector associated with element "iel"
%     iel - element number whose system dofs are to be determined
%     nnel - number of nodes per element
%     ndof - number of dofs per node
%-----------------------------------------------------------

edof = nnel*ndof;
k=0;
for i=1:nnel
    start = (nd(i)-1)*ndof;
    for j=1:ndof
        k=k+1;
        index(k)=start+j;
    end
end
%
%--------------------------------------------------------------------
%
function [kk,ff]=feaplyc2(kk,ff,bcdof,bcval)
%  Purpose:
%     Apply constraints to matrix equation [kk]{x}={ff}
%
%  Synopsis:
%     [kk,ff]=feaplybc(kk,ff,bcdof,bcval)
%
%  Variable Description:
%     kk - system matrix before applying constraints
%     ff - system vector before applying constraints
%     bcdof - a vector containging constrained d.o.f
%     bcval - a vector containing contained value
%
%     For example, there are constraints at d.o.f=2 and 10
%     and their constrained values are 0.0 and 2.5,
%     respectively.  Then, bcdof(1)=2 and bcdof(2)=10; and
%     bcval(1)=1.0 and bcval(2)=2.5.
%-----------------------------------------------------------
n=length(bcdof);
sdof=size(kk);
for i=1:n
    c=bcdof(i);
    for j=1:sdof
        kk(c,j)=0;
    end 
    kk(c,c)=1;
    ff(c)=bcval(i);
end
%
%--------------------------------------------------------------------
%
function plot_truss(joint_def,member_def)
%
jsiz=size(joint_def);
msiz=size(member_def);
nj=jsiz(1);
nm=msiz(1);
x=joint_def(:,1);
y=joint_def(:,2);
ei=member_def(:,1);
ej=member_def(:,2);
icolor=1;
%
%
A=member_def(:,3);
Amax=max(A);
%
labels=1;
if labels==1
    if icolor==1
        for i=1:nj
            plot(x(i),y(i),'*') % Drawing nodes
            t=[' ',num2str(i)];
            text(x(i)+0.015,y(i)+0.015,t,'Color','b');
            hold on
        end
    end
end
for k=1:nm
    if A(k)>.01*Amax
        str2=num2str(k);
        Xi(k)=x(ei(k));
        Yi(k)=y(ei(k));
        Xj(k)=x(ej(k));
        Yj(k)=y(ej(k));
        if icolor==1
            line([Xi(k) Xj(k)],[Yi(k) Yj(k)],'LineWidth',2,'Color','b')
        end
        if icolor==2
            line([Xi(k) Xj(k)],[Yi(k) Yj(k)],'LineWidth',2,'Color','r')
        end
        if icolor==3
            line([Xi(k) Xj(k)],[Yi(k) Yj(k)],'LineWidth',2,'Color','g')
        end
        if icolor==1
            if labels==1
                a=(Xi(k)+Xj(k))/2;
                b=(Yi(k)+Yj(k))/2;
                tt=[' e',num2str(k)];
                text(a,b,tt,'Color','r');
            end
        end
        hold off
    end
end
axis equal

