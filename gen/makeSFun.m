function makeSFun(model_path)

maxima_path= getenv('maxima_path');
if isempty(maxima_path)
    error('Please set path to maxima program via "setenv(''maxima_path'', ''PATH_TO_MAXIMA.EXE'');');
end
cagem_path= getenv('cagem_path');
if isempty(cagem_path)
    error('Please set path to cagem.mac script via "setenv(''cagem_path'', ''PATH_TO_CAMGEM.MAC'');');
end

cagem_base= fileparts(cagem_path);
cagem_base= fileparts(cagem_base);

[~, model_name]= fileparts(model_path);
target_path= fullfile(pwd, model_name);

if ~exist(target_path, 'dir')
    mkdir(target_path);
end

skip_gen= false;
code_file= fullfile(target_path, [model_name 'System2.hpp']);
if exist(code_file, 'file')
    dd_src= dir(model_path);
    dd_dst= dir(code_file);
    if dd_src.datenum<dd_dst.datenum
        skip_gen= true;
    end
end

if ~skip_gen
    command_str= [
        maxima_path, ...
        ' ', ...
        '--batch-string="', ...
        'load(\"' cagem_path '\")\$ ', ...
        'sys:cagem(\"' model_path '\")\$ ', ...
        'cagem2c2(sys, \"' target_path '\")\$', ...
        '"'];
        
    fprintf('%s\n', command_str);
        
    [~, res]= system(command_str);
    fprintf('Maxima said:\n%s\n', res);
end

generated= false;
if exist(code_file, 'file')
    dd_src= dir(model_path);
    dd_dst= dir(code_file);
    if dd_src.datenum<dd_dst.datenum
        generated= true;
    end
end
if ~generated
    error('Code was not properly generated by Maxima');
end

start_dir= pwd;
cleanupVar= onCleanup(@() cd(start_dir));

cd(target_path);

sfun_name= [model_name '_sfun'];

clear mex
output_str= '-output';

sfun_file= fullfile(cagem_base, 'gen', 'CADyn_sfun.cpp');

defines= {['MBSystem=' model_name]};
defines{end+1}= ['S_FUNCTION_NAME=' sfun_name];

defines= strcat('-D', defines);

options= {'-g' output_str sfun_name};
options= [options {'CXXFLAGS="$CXXFLAGS -std=c++11 -Wall -fdiagnostics-show-option"'}];
    
includes= {
    '.'
    fullfile(cagem_base, 'src')
    };

includes= strcat('-I', includes);

sources= {
    fullfile(cagem_base, 'src', 'ODEOrder2.cpp')
    };

sfun_name_ext= [sfun_name '.' mexext];
if exist(fullfile(start_dir, sfun_name_ext), 'file')
    delete(fullfile(start_dir, sfun_name_ext));
end

mex(options{:}, defines{:}, includes{:}, sfun_file, sources{:});

compiled= false;
if exist(sfun_name_ext, 'file')
    dd_src= dir(code_file);
    dd_dst= dir(sfun_name_ext);
    if dd_src.datenum<dd_dst.datenum
        compiled= true;
    end
end
if ~compiled
    error('sfun was not properly compiled');
end

movefile(sfun_name_ext, start_dir);