old_dir= pwd;
cd('../../../gen');
gen_dir= pwd;
addpath(gen_dir);
cd(old_dir)
setenv('maxima_path', '/usr/bin/maxima')
setenv('cagem_path', fullfile(gen_dir, 'cagem.mac'))
makeMex('../../turbine_coll_flap_edge_pitch.mac')
