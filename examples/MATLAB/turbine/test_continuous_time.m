x0= [0.1 0.01 0 1 0.01 0]';
dx0= [0 0 1 0 0 1.001]';
u= [0 0 0 0 0 0]';

[x1, dx1, ~, S]= turbine_coll_flap_edge_pitch_mex(x0, dx0, u, p);

E= blkdiag(eye(6), S(1:6,1:6));
A= S(:, 1:12);
A(1:6, 1:6)= 0;
B= S(:, 13:18);
