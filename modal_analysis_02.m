% MODAL ANALYSIS FOR CONTINUOUS INPUT
clear all;

% call function
output_list = modal_analysis_2();

% define function
function output_list = modal_analysis_2()
    % define parameters
    num_masses = 300;
    total_mass = 2;
    tension_force = 5;
    string_length = 5000;
    damping_coeff = 0;
    dx = string_length/(num_masses+1);
    w = 2;
    h = 10;
    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    Uf_func = @(t) wave_func_in(t,w,h);
    dUfdt_func = @(t) wave_func_derivative(t,w,h);
    %generate the struct
    string_params = struct();
    string_params.n = num_masses;
    string_params.M = total_mass;
    string_params.Uf_func = Uf_func;
    string_params.dUfdt_func = dUfdt_func;
    string_params.Tf = tension_force;
    string_params.L = string_length;
    string_params.c = damping_coeff;
    string_params.dx = dx;
    %initial conditions
    U0 = zeros(num_masses,1);
    dUdt0 = zeros(num_masses,1);
    V0 = [U0;dUdt0];
    tspan = [0,75];

    c = sqrt(string_params.Tf/(total_mass/string_params.L));
    num_frequencies = 10;
    Bn = 5;
    n = transpose(1:num_frequencies);
    resonant_frequencies = c*pi*n/string_params.L;
    mode_shapes = zeros(num_frequencies,1);
    for x = 1:(string_params.L-1)/10:string_params.L
        mode_shapes(:,end+1) = Bn*sin(pi*n*x/string_params.L)
    end

    output_list = [mode_shapes,resonant_frequencies];
end


