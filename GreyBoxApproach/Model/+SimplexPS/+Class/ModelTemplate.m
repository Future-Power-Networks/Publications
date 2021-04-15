% This class shows a template for models based on "ModelAdvance". 

% Author(s): Yitong Li

%% Notes
%
% This is just a template. For practical examples, please see
% "Inductor.m", which is a single-phase inductor, and see
% "SynchronousMachine.m", which is a three-phase synchronous machine.

%% Class

% Please change "ModelTemplate" here to your customized name and save this
% file as your customized name as well.
classdef ModelTemplate < SimplexPS.Class.ModelAdvance
    
  	methods
        % Constructor
        % Please change "ModelTemplate" to your customized name.
        function obj = ModelTemplate(varargin)
            setProperties(obj,nargin,varargin{:});
        end
    end
    
    methods(Static)
        
        % Set the strings of state, input, output
        % These strings are mainly for printing output and facilitate the
        % users.
        function [State,Input,Output] = SignalList(obj)
        	State  = {'x1','x2'}; 	% x, state
            Input  = {'v'};        	% u, input
            Output = {'i'};        	% y, output
        end
        
        % Calculate the equilibrium
        % The equilibrium is determined by the power flow data and device's
        % own paramters.
        function [x_e,u_e,xi] = Equilibrium(obj)
         	% Get the power PowerFlow values
            P 	= obj.PowerFlow(1);
            Q	= obj.PowerFlow(2);
            V	= obj.PowerFlow(3);
            xi	= obj.PowerFlow(4);
            w   = obj.PowerFlow(5);
            
            % Get parameters
            obj.Para(1);
            
            % Calculate equilibrium
            
            % Set equilibrium
            x_e = [];
            u_e = [];
            xi  = xi;
        end
        
    	% State space model
        % This function defines the state space model of this device,
        % and is the core part for capturing the dynamics of this device.
        %
        % The state space model should be a large-signal model rather than
        % a small-signal model. The linearized model will be calculated by
        % functions in the parent class and the linearization point (i.e.
        % equilibrium) is calculated above.
        function [Output] = StateSpaceEqu(obj,x,u,CallFlag)
          	% Get parameter
            obj.Para(1);
            
        	% Get state
            x(1);
            
            % Get input
            u(1);
            
            % State space equations
          	% dx/dt = f(x,u)
            % y     = g(x,u)
            if CallFlag == 1
                % ### Call state equation: dx/dt = f(x,u)
                
                f_xu = [];
                Output = f_xu;
            elseif CallFlag == 2
                % ### Call output equation: y = g(x,u)
                
                g_xu = [];
                Output = g_xu;
            end
        end
        
    end
end