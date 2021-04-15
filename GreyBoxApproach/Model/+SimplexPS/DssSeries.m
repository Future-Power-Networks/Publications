% This function achieves the "series" operation of two descriptor state
% space systems.

% Author(s): Yitong Li

%% Notes

% Schematic:
% u1 -> Gdss1 -> y1 = u2 -> Gdss2 -> y2

%% Example

% Gdss1:
% E1*dx1/dt = A1*x1 + B1*u1
% y1        = C1*x1 + D1*u1
% Gdss2:
% E2*dx2/dt = A2*x2 + B2*u2
% y2        = C2*x2 + D2*u2
% ->
% Gdss
% [E2   ]*[dx2]/dt = [A]*[x2] + [B]*[u1]
% [   E1] [dx1]          [x1]
% y2               = [C]*[x2] + [D]*[u1]
%                        [x1]
% In other words, the series operation does not influence the E matrix,
% which is valid at least when ly1 = lu2.

%%
function Gdss = DssSeries(Gdss1,Gdss2)

    if (~is_dss(Gdss1)) || (~is_dss(Gdss2))
        error(['Error: System 1 and/or 2 is not in dss form.']);
    end
    [lx1,lu1,ly1] = dss_GetDim(Gdss1);
    [lx2,lu2,ly2] = dss_GetDim(Gdss2);
    
    if ly1 ~= lu2
        error(['Error: The output dimention of system 1 does not equal to the the input dimension of system 2.']);
    end
    
    E1 = Gdss1.E;
    E2 = Gdss2.E;
    
    Gss1 = dss2ss(Gdss1);
    Gss2 = dss2ss(Gdss2);
    
    Gss = series(Gss1,Gss2);
    Gdss = Gss;
    Gdss.E = blkdiag(E2,E1);
    
end