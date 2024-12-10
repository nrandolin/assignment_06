% MODAL ANALYSIS FOR CONTINUOUS INPUT
clear all;

% call function
output_list = modal_analysis_2();

mode_num = 1;
mode_shape = output_list(:,1+mode_num);
resonant_freq = output_list(mode_num, end-1);
frequencies = output_list(mode_num, end);



%% define function
function output_list = modal_analysis_2()
    % define parameters
    num_masses = 15;
    total_mass = 30;
    tension_force = 8;
    string_length = 3;
    damping_coeff = 0.1;
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
    tspan = [0,1];
    x_pos = dx.*[1:1:num_masses]
    x_pos = [0, x_pos, 3]

    c = sqrt(string_params.Tf/(total_mass/string_params.L));
    num_frequencies = 10;
    Bn = 5;
    n = transpose(1:num_frequencies);
    frequencies = pi*n/string_params.L;
    resonant_frequencies = c*pi*n/string_params.L;
    mode_shapes = ones(num_frequencies,1);
    for i = 1:length(x_pos)
        mode_shapes(:,end+1) = Bn*sin(pi*n.*x_pos(i)/string_params.L);
    end
    amplitude_Uf = 5;
    mode_shapes = mode_shapes(:,2:end)
    output_list = [mode_shapes,resonant_frequencies, frequencies];
    %string_simulation_02(num_masses, total_mass, tension_force,...
    %string_length, damping_coeff, dx,amplitude_Uf,frequencies(2,1))
    figure();
    plot(n, resonant_frequencies, '.b', 'MarkerSize', 8)
    title("Resonant Frequency vs. n")
    hold off
    figure();
    plot(x_pos, mode_shapes(2,:), '.r', 'MarkerSize', 20)
    hold on
    plot(x_pos, mode_shapes(2,:), '-b')
end

%%
%construct the nxn discrete laplacian matrix
% n = 5;
% I_n = eye(n); % build the nxn identity matrix
% my_Laplacian = circshift(I_n, [0,1]) - 2*I_n + circshift(I_n, [0,-1])
% my_Laplacian(1,n) = 0; 
% my_Laplacian(n,1) = 0;
%Use MATLAB to solve the generalized eigenvalue problem
function string_simulation_02(num_masses, total_mass, tension_force,...
    string_length, damping_coeff, dx,amplitude_Uf,omega_Uf)
    %list of x points (including the two endpoints)
    xlist = linspace(0,string_length,num_masses+2);
    Uf_func = @(t_in) amplitude_Uf*cos(omega_Uf*t_in);
    dUfdt_func = @(t_in) -omega_Uf*amplitude_Uf*sin(omega_Uf*t_in);
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
    tspan = [0,10];
    %run the integration
    ralston_struct.A = [0, 0, 0; 0.5, 0, 0; 0, 0.75, 0];
    ralston_struct.B = [2/9, 1/3, 4/9];
    ralston_struct.C = [0; 0.5; 0.75];
    h_ref = 0.1;
    my_rate_func = @(t,V) string_rate_func01(t,V,string_params);
    
    [tlist,Vlist, ~, ~] = explicit_RK_fixed_step_integration ...
    (my_rate_func,tspan,V0,h_ref,ralston_struct);
    
    figure(1);

    axis manual;
    for i = 1:length(tlist)
%         for mass_num = 1:length(Vlist(1,:))/2
            x_pos = dx.*[1:1:num_masses];
            position_of_masses = Vlist(i,1:num_masses);
            plot(x_pos,position_of_masses,"b.","MarkerSize",20);
            hold on;
            plot(max(x_pos)+dx,Uf_func(tlist(i)),"r.","MarkerSize",20);
            hold on;
            plot(0,0,"r.","MarkerSize",20);
            x_pos_all = [0, x_pos, max(x_pos)+dx];
            position_of_masses_all = [0, position_of_masses, Uf_func(tlist(i))];
            hold on;
            plot(x_pos_all,position_of_masses_all,"b-");
            ylim([-15,15]);
            xlabel("X")
            ylabel("Y")
            title("Mode Shape Wave Animation")
%         end
        pause(0.1);
        clf;
    end

    figure(2);
    axis manual;
    i = 100;
        x_pos = dx.*[1:1:num_masses];
        position_of_masses = Vlist(i,1:num_masses);
        plot(x_pos,position_of_masses,"b.","MarkerSize",20);
        hold on;
        plot(max(x_pos)+dx,Uf_func(tlist(i)),"r.","MarkerSize",20);
        hold on;
        plot(0,0,"r.","MarkerSize",20);
        x_pos_all = [0, x_pos, max(x_pos)+dx];
        position_of_masses_all = [0, position_of_masses, Uf_func(tlist(i))];
        hold on;
        plot(x_pos_all,position_of_masses_all,"b-");
        ylim([-15,15]);
        xlabel("X")
        ylabel("Y")
        title("Mode Shape")

end
%% explicit_RK fixed_step integrator
%Runs numerical integration arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%tspan: a two element vector [t_start,t_end] that denotes the integration endpoints
%X0: the vector describing the initial conditions, X(t_start)
%h_ref: the desired value of the average step size (not the actual value)
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%t_list: the vector of times, [t_start;t_1;t_2;...;.t_end] that X is approximated at
%X_list: the vector of X, [X0’;X1’;X2’;...;(X_end)’] at each time step
%h_avg: the average step size
%num_evals: total number of calls made to rate_func_in during the integration
function [t_list,X_list,h_avg, num_evals] = explicit_RK_fixed_step_integration ...
(rate_func_in,tspan,X0,h_ref,BT_struct)
    % calculate steps and h
    [num_steps, h_avg] = iteration_solver(tspan, h_ref);
    % define variables
    XA = X0;
    num_evals = 0;
    t_list = linspace(tspan(1),tspan(2),num_steps+1);
 
    X_list = zeros(num_steps+1,length(X0));
    X_list(1,:) = X0';
    %calculate the values until it is just short of the end value
    for i = 1:num_steps
        t = t_list(i);
        [XB, temp_eval] = explicit_RK_step(rate_func_in,t,XA,h_avg,BT_struct);
        num_evals = num_evals + temp_eval;

        X_list(i+1,:)= XB';
        XA = XB;
    end  
end

%% explicit_RK_step
%This function computes the value of X at the next time step
%for any arbitrary RK method
%INPUTS:
%rate_func_in: the function used to compute dXdt. rate_func_in will
% have the form: dXdt = rate_func_in(t,X) (t is before X)
%t: the value of time at the current step
%XA: the value of X(t)
%h: the time increment for a single step i.e. delta_t = t_{n+1} - t_{n}
%BT_struct: a struct that contains the Butcher tableau
% BT_struct.A: matrix of a_{ij} values
% BT_struct.B: vector of b_i values
% BT_struct.C: vector of c_i values
%OUTPUTS:
%XB: the approximate value for X(t+h) (the next step)
% formula depends on the integration method used
%num_evals: A count of the number of times that you called
% rate_func_in when computing the next step
function [XB, num_evals] = explicit_RK_step(rate_func_in,t,XA,h,BT_struct)
    k = zeros(length(XA),length(BT_struct.B));
    for i = 1:length(BT_struct.B)
        k(:,i) = rate_func_in(t+BT_struct.C(i)*h, XA+h*(k*BT_struct.A(i,:)'));
    end
    XB = XA + h*(k*BT_struct.B');
    num_evals = length(BT_struct.B);
end
%% ITERATION SOLVER
function [num_steps, h] = iteration_solver(tspan, h_ref)
    t_range = tspan(2)-tspan(1);
    num_steps = t_range/h_ref;%The number of steps is the range devided by goal h 
    num_steps = ceil(num_steps);%Round the number of steps up (to get a real number)
    h = t_range/num_steps; % Divide range by steps to get real h
end

