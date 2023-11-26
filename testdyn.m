clear
clc

alpha=[0,     pi/2,  0,           0,         pi/2,      -pi/2];
a=    [0,     0,     -300,      -276,    0,         0].*0.001;
d=    [121.5, 0,     0,           110.5,    90,     82].*0.001;
theta=[0,     -pi/2,  0,           -pi/2,     0,         0];
dh = [alpha; a; d; theta];

myRob=mandoRob(dh,'MDH');
rob=[0,0,0,0,0,0];%0 standard for rotation, 1 standard for translation
N = size(rob);
g = sym('g');
assume(g,'real');
m = sym('m',N); 
assume(m, 'real')
Pc = sym('Pc%d%d',[3,6]);
assume(Pc,'real');
Ic = sym('Ic_%d%d_%d',[3,3,6]);
assume(Ic,'real');
for i = 1:size(Ic,3)
    Ic(2,1,i)=Ic(1,2,i);
    Ic(3,1,i)=Ic(1,3,i);
    Ic(3,2,i)=Ic(2,3,i);
end
q = sym('q',N);       % 1st group of state variables (joints' angle or position)
assume(q, 'real')
dq = sym('dq',N);     % 2nd group of state variables (joints' velocity)
assume(dq, 'real')

mass_list = [2.920; 6.787; 2.450; 1.707; 1.707; 0.176]; 
%连杆质心相对于坐标系{i}坐标(m)
mass_center_list = [0.0316 -3.1464 -13.8983;
                    131.5620 -0.0210 112.1840;
                    190.3840 0.0410 17.1800;
                    0.0886 21.0083 -2.5014;
                    -0.0886 -21.0083 -2.5014;
                    0 0 8.0000]*10^-3;
mass_center_list = mass_center_list';
%inertia_tensor_list 连杆相对于质心坐标系的惯量张量(kg*m^2)
inertia_tensor_list = zeros(3,3,6);
inertia_tensor_list(:,:,1)  = [42.614 0.046 0.062; 0.046 41.164 -1.386; 0.062 -1.386 31.883]*10^-4;
inertia_tensor_list(:,:,2)  = [100.7 -1.8 1.6; -1.8 1100.8 0; 1.6 0 1087.1]*10^-4;
inertia_tensor_list(:,:,3)  = [31.45 0.48 7.23; 0.48 172.41 -0.15; 7.23 -0.15 166.82]*10^-4;
inertia_tensor_list(:,:,4)  = [20.92 -0.061 0.078; -0.061 16.808 0.992; 0.078 0.992 19.75]*10^-4;
inertia_tensor_list(:,:,5)  = [20.92 -0.061 -0.078; -0.061 16.808 -0.992; -0.078 -0.992 19.75]*10^-4;
inertia_tensor_list(:,:,6)  = [0.9296 0 0; 0 0.9485 0; 0 0 1.5925]*10^-4; 
%% prepare substitute data 
sym_sub = [];
num_sub = [];
for i = 1:size(Ic,3)
    sym_sub=[sym_sub,m(i)];
    num_sub=[num_sub,mass_list(i)];
    for j = 1:size(Pc,1)
        sym_sub=[sym_sub,Pc(j,i)];
        num_sub=[num_sub,mass_center_list(j,i)];
    end
    for j = 1:size(Ic,1)
        for k = 1:size(Ic,2)
            sym_sub=[sym_sub,Ic(j,k,i)];
            num_sub=[num_sub,inertia_tensor_list(j,k,i)];
        end
    end
end
[J_L,M_L,C_L,G_L] = Lagrange_InverseDyn(dh,rob,Pc,Ic,g,m,q,dq);
G = subs(G_L,[g,sym_sub],[9.81,num_sub]);
M = subs(M_L,sym_sub,num_sub);
C = subs(C_L,sym_sub,num_sub);



[MM,CC,GG]=myRob.lagrange_dyn(num_sub);

%% calculate tau
thetalist = [0;pi/3;pi/3;-pi/6;pi/2;pi/2];
dthetalist = [pi/2;pi/2;pi/6;0;pi/2;pi/4];
ddthetalist = [0;pi/34;0;pi/5;0;0.34];

tau_Rob=(double(subs(MM,q,(thetalist)'))*ddthetalist+double(subs(CC,[q,dq],[(thetalist)',dthetalist']))*dthetalist + double(subs(GG,q,thetalist')))';
tau_L = (double(subs(M,q,(thetalist)'))*ddthetalist+double(subs(C,[q,dq],[(thetalist)',dthetalist']))*dthetalist + double(subs(G,q,(thetalist'))))';



