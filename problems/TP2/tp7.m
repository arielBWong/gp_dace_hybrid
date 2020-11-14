classdef tp7
    properties
        p = 1;
        q = 1;
        n_lvar;
        n_uvar;
        xu_bl;
        xu_bu;
        xl_bl;
        xl_bu;
        xprime;
        name = 'tpso7';
        uopt = NaN;
        lopt = 0.0625; % double check needed
    end
    methods
        function obj = tp7(p, q)
            % level variables
             if nargin > 1
                 obj.q = q;
                 obj.p = p;
             end
            obj.n_lvar = obj.q;
            obj.n_uvar = obj.p;
            obj.xprime = ones(1,  obj.q)  * 8;
            % bounds
            %init bound upper level
            obj.xu_bl = zeros(1, obj.p);
            obj.xu_bu = ones(1, obj.p);
           
           
            % init bound lower level
            obj.xl_bl = zeros(1, obj.q) ;
            obj.xl_bu = ones(1, obj.q) * 10;      
        end
        
        function [f, c] = evaluate_u(obj, xu, xl) 
            %-obj
            f = [];
            c = [];
            
        end
        
        
        function [f, c] = evaluate_l(obj, xu, xl)
           
            %-obj
            f = 0;
            n = 1;
           for i = 1: obj. q
               index  = xl(:, i) > 2 & xl(:, i) <= 8;
               other = ~index;
               
               fb1  = (xl(:, i) - 2).^4 ./ 6;
               scale = (8-2)^4/6;
               fb1 = fb1 ./ scale;
               
               fb1(other) = 0;
               f = f+ fb1 ;
           end
           
           f = f ./ n;
           f = - f;
            %-con
            c = [];
         end
    end
end
