
clc;clear
alpha=[0,     pi/2,  0,           0,         pi/2,      -pi/2];
a=    [0,     0,     -300,      -276,    0,         0].*0.001;
d=    [121.5, 0,     0,           110.5,    90,     82].*0.001;
theta=[0,     -pi/2,  0,           -pi/2,     0,         0];
dh = [alpha; a; d; theta];

% thetalist = [0;pi/3;pi/3;-pi/6;pi/2;pi/2];
thetalist = [0;pi/3;pi/3;-pi/3;pi/3;pi/3];

myRob=mandoRob(dh,'MDH');
an=myRob.Rotate(thetalist')
an=double(an);


alph = dh(1,:);
a = dh(2,:);
d = dh(3,:);
theta = dh(4,:);
L1 = Revolute('d', d(1), 'a', a(1), 'alpha', alph(1), 'modified');
L2 = Revolute('d', d(2), 'a', a(2), 'alpha', alph(2), 'modified');
L3 = Revolute('d', d(3), 'a', a(3), 'alpha', alph(3), 'modified');
L4 = Revolute('d', d(4), 'a', a(4), 'alpha', alph(4), 'modified');
L5 = Revolute('d', d(5), 'a', a(5), 'alpha', alph(5), 'modified');
L6 = Revolute('d', d(6), 'a', a(6), 'alpha', alph(6), 'modified');

robot=SerialLink([L1,L2,L3,L4,L5,L6],'name','6dof_arm','comment','LL');  %SerialLink类函数

ann=robot.fkine(theta+thetalist');

