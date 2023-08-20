function code_file= target_file_name(model_name, file_to_generate)
switch file_to_generate
    case 'cpp_direct'
        code_file= [model_name '_direct.hpp'];
    case 'cpp_direct_gmres'
        code_file= [model_name '_direct_gmres.hpp'];
    case 'acados'
        code_file= [model_name '_acados.m'];
    otherwise
        error('Unknown file to generate "%s"', file_to_generate)
end
