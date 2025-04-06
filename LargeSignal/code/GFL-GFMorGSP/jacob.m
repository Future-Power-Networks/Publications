function [A,B,C,D] = jacob(f,g,x_e,u_e)

    % Calculate equilibrium of dx_e and y_e
    % The function "StateSpaceEqu" is defined in the subclass
    dx_e = f(x_e, u_e);
    y_e  = g(x_e, u_e);

    % Calculate length
    lx = length(x_e);
    lu = length(u_e);
    ly = length(y_e);

    % Initialize A,B,C,D
    A = zeros(lx,lx);
    B = zeros(lx,lu);
    C = zeros(ly,lx);
    D = zeros(ly,lu);

    % Get the perturb size
    perturb_factor = 1e-6;

    % Perturb x to calculate Ass and Css
    for i = 1:lx
        x_p = x_e;                      % Reset xp
        perturb_x = perturb_factor * abs(1+abs(x_e(i))); 	% Perturb size
        x_p(i) = x_e(i) + perturb_x;                        % Positive perturb on the ith element of xp
        dx_p = f(x_p, u_e);
        y_p  = g(x_p, u_e);
        A(:,i) = (dx_p - dx_e)/(x_p(i) - x_e(i));
        C(:,i) = (y_p - y_e)/(x_p(i) - x_e(i));
    end

    % Perturb u to calculate Bss and Dss
    for i = 1:lu
        up = u_e;                       % Reset up
        perturb_u = perturb_factor * abs(1+abs(u_e(i)));    % Perturb size
        up(i) = u_e(i) + perturb_u;                         % Positve perturb on the ith element of up
        dx_p = f(x_p, u_e);
        y_p  = g(x_p, u_e);
        B(:,i) = (dx_p - dx_e)/(up(i) - u_e(i));
        D(:,i) = (y_p - y_e)/(up(i) - u_e(i));
    end 
end  