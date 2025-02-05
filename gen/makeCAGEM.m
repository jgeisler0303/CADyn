function makeCAGEM(model_path, target_path, skip_gen, files_to_generate, reduce_consts)
maxima_path= getenv('maxima_path');
if isempty(maxima_path)
    error('Please set path to maxima program via "setenv(''maxima_path'', ''PATH_TO_MAXIMA.EXE'');');
end

cagem_path= getenv('cagem_path');
if isempty(cagem_path)
    error('Please set path to cagem.mac script via "setenv(''cagem_path'', ''PATH_TO_CAMGEM.MAC'');');
end

if ~iscell(files_to_generate)
    files_to_generate= {files_to_generate};
end
if ~exist('reduce_consts', 'var')
    reduce_consts= 0;
end
if length(reduce_consts)>1
    const_matrix_once= reduce_consts(2);
    reduce_consts= reduce_consts(1);
else
    const_matrix_once= 0;
end
if reduce_consts
    warning('Reduce constants seems to be broken. It will therefore not be applied.')
    reduce_consts= 0;
end

[~, model_name]= fileparts(model_path);
if ~exist('target_path', 'var')
    target_path= fullfile(pwd, model_name);
end

if ~exist(target_path, 'dir')
    mkdir(target_path);
end

skip_gen_file= false(length(files_to_generate), 1);
if ~exist('skip_gen', 'var') || isempty(skip_gen)
    dd_src= dir(model_path);
    for i= 1:length(files_to_generate)
        code_file= get_target_name(model_name, files_to_generate{i});
        code_path= fullfile(target_path, code_file);
        if exist(code_path, 'file')
            dd_dst= dir(code_path);
            if dd_src.datenum<dd_dst.datenum
                fprintf('Skipping generation of "%s"\n', code_file)
                skip_gen_file(i)= true;
            end
        end
    end
    skip_gen= all(skip_gen_file);
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
        'sys:cagem(\"' model_path '\", ' num2str(reduce_consts) ')\$ ', ...
        ];

    template_gen_loaded= false;
    for i= 1:length(files_to_generate)
        if skip_gen_file(i), continue, end

        [gen_cmd, load_cmd]= get_generator(files_to_generate{i}, target_path, const_matrix_once);
        if isempty(gen_cmd), continue, end

        if isempty(load_cmd) && ~template_gen_loaded
            template_gen_loaded= true;
            command_str= [command_str, 'load(\"' strrep(getenv('cagem_path'), 'cagem.mac', 'cadyn_gen_template.mac') '\")\$ '];
        else
            command_str= [command_str, load_cmd];
        end
        
        command_str= [command_str, gen_cmd];
    end
    command_str= [command_str '"'];

    fprintf('%s\n', command_str);
        
    [status, res]= system(command_str);
    fprintf('Maxima said:\n%s\n', res);

    for i= 1:length(files_to_generate)
        if skip_gen_file(i), continue, end
        [final_name, code_file]= get_target_name(model_name, files_to_generate{i});
        code_path= fullfile(target_path, code_file);
        generated= false;
        if exist(code_path, 'file')
            dd_src= dir(model_path);
            dd_dst= dir(code_path);
            if dd_src.datenum<dd_dst.datenum
                generated= true;
            end
        end
        if ~generated
            error('Code file "%s" was not properly generated by Maxima', code_file);
        else
            if ~strcmp(code_file, final_name)
                movefile(fullfile(target_path, code_file), fullfile(target_path, final_name))
            end
        end
    end
else
    fprintf('Skipping CADyn generation because model is already up to date or you wanted it so\n');
end


function [final_name, target_name]= get_target_name(model_name, generator)
if generator(1)=='_'
    target_name= [model_name generator];
    final_name= target_name;
else
    target_name= [model_name '_' generator];
    final_name= generator;
end

function [gen_cmd, load_cmd]= get_generator(generator, target_path, const_matrix_once)
if generator(1)=='_', generator(1)= []; end

template_name= ['cadyn_' generator '.tem'];
if exist(strrep(getenv('cagem_path'), 'cagem.mac', fullfile('templates', template_name)), 'file')
    load_cmd= '';
    gen_cmd= ['cadyn_gen_template(sys, \"' template_name '\", \"' target_path '\", ' num2str(const_matrix_once) ')\$ '];
else
    gen_script= ['cadyn_gen_' generator];
    gen_script_file= [gen_script '.mac'];
    [~, gen_script]= fileparts(gen_script);
    if ~exist(strrep(getenv('cagem_path'), 'cagem.mac', gen_script_file), 'file')
        warning('Don''t know how to make %s. No template and no script found', generator)
        load_cmd= '';
        gen_cmd= '';
    else
        load_cmd= ['load(\"' strrep(getenv('cagem_path'), 'cagem.mac', gen_script_file) '\")\$ '];
        gen_cmd= [gen_script '(sys, \"' target_path '\")\$ '];
    end
end
