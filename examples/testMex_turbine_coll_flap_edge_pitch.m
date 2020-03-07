old_dir= pwd;
cd('../gen');
gen_dir= pwd;
addpath(gen_dir);
cd(old_dir)
setenv('maxima_path', '/usr/bin/maxima')
setenv('cagem_path', fullfile(gen_dir, 'cagem.mac'))
makeMex('turbine_coll_flap_edge_pitch.mac')

%%
p.BlCrEdg= -50000;
p.BlCrFlp= -50000;
p.BlCtEdg= 1000;
p.BlCtFlp= 1000;
p.BlDEdg= 60;
p.BlDFlp= 20;
p.BlIner= 20000000;
p.BlKEdg= 30000;
p.BlKFlp= 7000;
p.BlMass= 10000;
p.BlMeEdg= 700;
p.BlMeFlp= 500;
p.DTTorDmp= 2000000;
p.DTTorSpr= 600000000;
p.GBRatio= 100;
p.GenIner= 300;
p.HubIner= 30000;
p.HubMass= 30000;
p.NacMass= 100000;
p.TwD= 4000;
p.TwK= 600000;
p.TwMe= 30000;
p.TwTrans2Roll= 0.02;

x0= [0.1 0.01 0 1 0.01 0]';
dx0= [0 0 1 0 0 1.001]';
u= [0 0 0 0 0 0];

ts= 0.01;
t= 0:ts:10;
nt= length(t);

%%
x= zeros(length(x0), nt);
dx= zeros(length(x0), nt);
x(:, 1)= x0;
dx(:, 1)= dx0;
cpu_time= zeros(1, nt);
int_err= zeros(1, nt);
n_steps= zeros(1, nt);
n_backsteps= zeros(1, nt);

for i= 2:nt
    [x(:, i), dx(:, i), ~, ~, cpu_time(i), int_err(i), n_steps(i), n_backsteps(i)]= turbine_coll_flap_edge_pitch_mex(x(:, i-1), dx(:, i-1), u, p, ts, struct('StepTol', 1e-1, 'AbsTol', 1e-4, 'RelTol', 1e-4));
end


%%
clf
subplot(3, 1, 1)
plot(t, cpu_time)
grid on
xlabel('time in s')

subplot(3, 1, 2)
plot(t, int_err)
grid on
xlabel('time in s')

subplot(3, 1, 3)
plot(t, n_steps, t, n_backsteps)
grid on
xlabel('time in s')
