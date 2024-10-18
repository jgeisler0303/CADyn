function code_file= target_file_name(model_name, file_to_generate)
switch file_to_generate
    case 'matlab'
        code_file= [model_name '_matlab_nonlin.m'];
    case 'ocp_direct'
        code_file= ['ipddp_problem_' model_name '.hpp'];
    case 'ocp_direct_rk1i'
        code_file= ['ipddp_problem_' model_name '_rk1i.hpp'];
    case 'ocp_direct_rk1'
        code_file= ['ipddp_problem_' model_name '_rk1.hpp'];
    case 'descriptor_form'
        code_file= [model_name '_descriptor_form.hpp'];
    case 'cpp_direct'
        code_file= [model_name '_direct.hpp'];
    case 'hpp_direct_rk1i'
        code_file= [model_name '_direct_rk1i.hpp'];
    case 'cpp_direct_gmres'
        code_file= [model_name '_direct_gmres.hpp'];
    case 'acados'
        code_file= [model_name '_acados.m'];
    otherwise
        error('Unknown file to generate "%s"', file_to_generate)
end
