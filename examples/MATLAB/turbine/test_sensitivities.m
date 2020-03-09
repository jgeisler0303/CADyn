x0= [0.1 0.01 0 1 0.01 0]';
dx0= [0 0 1 0 0 1.001]';
u= [0 0 0 0 0 0]';

ts= 0.01;
options= struct('StepTol', 1e-1, 'AbsTol', 1e-4, 'RelTol', 1e-4);

% ts= 0.1;
% options= struct('StepTol', 1e1, 'AbsTol', 1e-2, 'RelTol', 1e-2);

[x1, dx1, ~, S]= turbine_coll_flap_edge_pitch_mex(x0, dx0, u, p, ts, options);

xdxu= [x0; dx0; u];

S_= zeros(12, 18);
for i= 1:18
    xdxu_= xdxu;
    xdxu_(i)= xdxu_(i)+1e-5;
    [x1_, dx1_]= turbine_coll_flap_edge_pitch_mex(xdxu_(1:6), xdxu_(7:12), xdxu_(13:18), p, ts, options);
    S_(:, i)= ([x1_; dx1_] - [x1; dx1])/1e-5;
end

norm(S-S_)

