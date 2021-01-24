function makeMex(model_name, target_path)

cagem_path= getenv('cagem_path');
if isempty(cagem_path)
    error('Please set path to cagem.mac script via "setenv(''cagem_path'', ''PATH_TO_CAMGEM.MAC'');');
end

cagem_base= fileparts(cagem_path);
cagem_base= fileparts(cagem_base);

start_dir= pwd;
cleanupObj = onCleanup(@()cd(start_dir));

cd(target_path);

mex_name= [model_name '_mex'];

v= ver;
is_matlab= ~strcmp(v(1).Name, 'Octave');
if is_matlab
    clear mex
    output_str= '-output';
else
    clear(mex_name);
    output_str= '-o';
end


mex_file= fullfile(cagem_base, 'gen', 'CADyn_mex.cpp');

defines= {['MBSystem=' model_name]};
%  if exist('sfun', 'var') && sfun
%      defines{end+1}= ['S_FUNCTION_NAME=' mex_name];
%  end

defines= strcat('-D', defines);

options= {'-g' output_str mex_name};
if is_matlab
    options= [options {'CXXFLAGS="$CXXFLAGS -std=c++11 -Wall -fdiagnostics-show-option"'}];
end
    
includes= {
    '.'
    fullfile(cagem_base, 'src')
    };

includes= strcat('-I', includes);

sources= {
    fullfile(cagem_base, 'src', 'ODEOrder2.cpp')
    };

%  libs= { % order matters!!!
%      };

%  lib_dirs= strcat('-L', {acados_lib_dir});
%  libs= strcat('-l', libs);

mex_name_ext= [mex_name '.' mexext];
if exist(fullfile(start_dir, mex_name_ext), 'file')
    delete(fullfile(start_dir, mex_name_ext));
end

if ~is_matlab
    old_cxxflags= getenv('CXXFLAGS');
    [status, cxxflags]= system('mkoctfile --print  CXXFLAGS');
    cxxflags= [cxxflags ' -std=c++11 -Wall -fdiagnostics-show-option '];
    cxxflags= strrep(cxxflags, char(10), ' ');
    setenv('CXXFLAGS', cxxflags);
end

mex(options{:}, defines{:}, includes{:}, mex_file, sources{:}); %, lib_dirs{:}, libs{:})

if ~is_matlab
    setenv('CXXFLAGS', old_cxxflags);
end

compiled= false;
if exist(mex_name_ext, 'file')
    dd_src= dir(code_file);
    dd_dst= dir(mex_name_ext);
    if dd_src.datenum<dd_dst.datenum
        compiled= true;
    end
end
if ~compiled
    error('Mex was not properly compiled');
end

movefile(mex_name_ext, start_dir);
