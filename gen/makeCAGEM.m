function makeCAGEM(model_path, target_path, skip_gen, files_to_generate)
maxima_path= getenv('maxima_path');
if isempty(maxima_path)
    error('Please set path to maxima program via "setenv(''maxima_path'', ''PATH_TO_MAXIMA.EXE'');');
end

cagem_path= getenv('cagem_path');
if isempty(cagem_path)
    error('Please set path to cagem.mac script via "setenv(''cagem_path'', ''PATH_TO_CAMGEM.MAC'');');
end

if ~exist('files_to_generate', 'var')
    files_to_generate= {'cpp_direct'};
end
if ~iscell(files_to_generate)
    files_to_generate= {files_to_generate};
end

cagem_base= fileparts(cagem_path);
cagem_base= fileparts(cagem_base);

[~, model_name]= fileparts(model_path);
if ~exist('target_path', 'var')
    target_path= fullfile(pwd, model_name);
end

if ~exist(target_path, 'dir')
    mkdir(target_path);
end

% TODO: check for file to be generated according to files_to_generate
code_file= fullfile(target_path, [model_name '_direct.hpp']);
if ~exist('skip_gen', 'var') || isempty(skip_gen)
    skip_gen= false;
    if exist(code_file, 'file')
        dd_src= dir(model_path);
        dd_dst= dir(code_file);
        if dd_src.datenum<dd_dst.datenum
            skip_gen= true;
        end
    end
end

if ~skip_gen
    if ispc % even on Windows the separator needs to be /
        cagem_path= strrep(cagem_path, '\', '/');
        model_path= strrep(model_path, '\', '/');
        target_path= strrep(target_path, '\', '/');
        model_name= strrep(model_name, '\', '/');
    end

    command_str= [
        maxima_path, ...
        ' ', ...
        '--batch-string="', ...
        'load(\"' cagem_path '\")\$ ', ...
        'sys:cagem(\"' model_path '\")\$ ', ...
        'load(\"' strrep(cagem_path, 'cagem.mac', 'cagem_gen_matlab.mac') '\")\$ ', ...
        'writeMATLAB_nonlinear(sys,\"' target_path '/' model_name '_matlab_nonlin.m' '\")\$ ', ...
        'writeMATLAB_linear(sys,\"' target_path '/' model_name '_matlab_lin.m' '\")\$ ', ...
        'load(\"' strrep(cagem_path, 'cagem.mac', 'cagem_gen_indices.mac') '\")\$ ', ...
        'cagem_indices(sys,\"' target_path '\")\$ ', ...
        ];
    for i= 1:length(files_to_generate)
        command_str= [command_str, ...
        'load(\"' strrep(cagem_path, 'cagem.mac', ['cagem_gen_' files_to_generate{i} '.mac']) '\")\$ ', ...
        'cagem_gen_' files_to_generate{i} '(sys, \"' target_path '\")\$ ', ...
        ];
    end
    command_str= [command_str '"'];

    fprintf('%s\n', command_str);
        
    [status, res]= system(command_str);
    fprintf('Maxima said:\n%s\n', res);
else
    fprintf('Skipping CADyn generation because model is already up to date\n');
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
