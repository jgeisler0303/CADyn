sys= ss(c2d(tf(1, [1 2 1]), 0.1));
sys_= sys;
sys_.b= [sys_.b eye(2)];
sys_.d= [0 0 0];

Q= eye(2);
R= eye(1);

[KEST,L,P]= kalman(sys_, Q, R, 'delayed');

x= [0; 0];
Sigma= eye(2); %P;
u= ones(1000, 1);
y_meas= ones(1000, 1);
X= [];
for i= 1:10
    X(:, i)= x;
    y= sys.c*x;
    x= sys.a*x + sys.b*u(i);
    
    K= sys.a*Sigma*sys.c'*(sys.c*Sigma*sys.c' + R)^-1;
%     Sigma= sys.a*Sigma*sys.a' + Q - K*sys.c*Sigma*sys.a';
    Sigma__= sys.a*Sigma*sys.a' + Q;
%     K= Sigma__*sys.c'*(sys.c*Sigma__*sys.c' + R)^-1;
    d= y_meas(i)-y;
    x= x + K*d; % should only be output in next iteration but input should be from current time step
    Sigma= (eye(length(x))-K*sys.c) * Sigma__;
    Sigma= 0.5*(Sigma+Sigma');
end
