function [Jee,M,C,G] = Lagrange_InverseDyn(dh,rob,Pc,Ic,g,m,q,dq)
%@input: dh stanard for modified DH, dh_list = [alpha; a; d; theta]; DH(4X6).
%@input: int rob the 0 standard for rotation, the 1 standard for translation. but also length(rob) equal size(dh,2);  rob(1X6).
%@input:Pc sandard the mess center of link, example: Pc(1,1) standsrd for the center of x axis on first link, Pc(2,1)standsrd for the center of y axis on first link,Pc(3,1)standsrd for the center of z axis on first link  Pc(3X6).
%@input:Ic sandard the inertia tensor of link,example:Ic(:,:,1)standard for the first link inertia tensor, and Ic(:,:,1) is a symmetry matrix   Ic(3X3X6).
%@input:g:9.81.
%@input:m:the mess of link,  m(1X6).
%@input:q:joint angle for every link,  q(1X6).
%@input:dq:joint angle velocity for every link,  dq(1X6).

%@output:Pee: the end position for robot, respect standard for the position of x y z axis.   Pee(3X1).
%@output:Ree: the end rotation for robot, it's is a rotation traslation matrix. Ree(3X3).
%@output:Jee: jacobians for robot contain linear and angle part, Jee(6X6).
%@output:M:M(q),  M(6X6).
%@output:C:C(q,dq),   C(6X6).
%@output:G:G(q),   G(1X6).
%% parameters
alpha = dh(1,:);
a = dh(2,:);
d = dh(3,:);
theta = dh(4,:);
N = size(dh,2); 

for i=1:N
   pc{i}=Pc(:,i);
end


%% R
Tn=eye(4);
for i=1:N
Rot_theta_q=[cos(q(i)) -sin(q(i)) 0 0;sin(q(i)) cos(q(i)) 0 0;0 0 1 0;0 0 0 1];   
Rot_theta=[cos(theta(i)) -sin(theta(i)) 0 0;sin(theta(i)) cos(theta(i)) 0 0;0 0 1 0;0 0 0 1];
Tran_d=[1 0 0 0; 0 1 0 0; 0 0 1 d(i); 0 0 0 1];
Tran_a=[1 0 0 a(i); 0 1 0 0; 0 0 1 0; 0 0 0 1];
Rot_alpha=[1 0 0 0; 0 cos(alpha(i)) -sin(alpha(i)) 0; 0 sin(alpha(i)) cos(alpha(i)) 0; 0 0 0 1];
T(:,:,i)=Rot_alpha*Tran_a*Tran_d*Rot_theta*Rot_theta_q;                            
R{i}=T(1:3,1:3,i);
p{i}=T(1:3,4,i);
Tn=Tn*T(:,:,i);
TT{i}=Tn;
end

for i=1:N              
    m_P{i}=TT{i}(1:3,4)+TT{i}(1:3,1:3)*pc{i};
end 

%% R0{i} stand for 0--> ith framework rotation transition matrix
R0{1}=R{1};
for i=2:N
    R0{i}=R0{i-1}*R{i};
end
Ree = R0{N};

%% linear part of jacobians
for i=1:N
    Jv{i}=jacobian(m_P{i},q);
end

%% angular part of jacobians
if rob(1)==0
    Jo{1}=[R0{1}(:,3),zeros(3,N-1)];
elseif rob(1)==1
    Jo{1}=zeros(3,N);
end

for i=2:N
    if rob(i)==1
        Jo{i}=zeros(3,N);
    elseif rob(i)==0
        Jo{i}=Jo{i-1};
        Jo{i}(:,i)=R0{i}(:,3);
    end
end
Jee = [Jv{N};Jo{N}];

%% M
for i = 1:size(rob,2)
    In{i} = Ic(:,:,i);
end

for i=1:N
    Mass{i}=m(i)*Jv{i}.'*Jv{i}+Jo{i}.'*R0{i}*In{i}*R0{i}.'*Jo{i};
end

M=0;
for i=1:N
    M=Mass{i}+M;
end

%% C
for k=1:N
   for s=1:N
       c(1)=.5*((diff(M(k,s),q(1))+diff(M(k,1),q(s))-diff(M(1,s),q(k)))*dq(1));
      for i=2:N
       c(i)=.5*((diff(M(k,s),q(i))+diff(M(k,i),q(s))-diff(M(i,s),q(k)))*dq(i))+c(i-1);
      end
      C(k,s)=c(N);
   end
end

%% G
P(1)=m(1)*[0,0,g]*m_P{1};
for i=2:N
    P(i)=P(i-1)+m(i)*[0,0,g]*m_P{i};
end
P=P(N);
for i=1:N
    G(i,:)=diff(P,q(i));
end

%% dynamic equations
%  M(q)*ddq + C(q,dq)dq + G(q) = u
end