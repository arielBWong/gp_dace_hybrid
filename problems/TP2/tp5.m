classdef tp5
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
        name = 'tpso5';
        uopt = NaN;
        lopt = NaN; % double check needed
    end
    methods
        function obj = tp5(p, q)
            % level variables
             if nargin > 1
                 obj.q = q;
                 obj.p = p;
             end
            obj.n_lvar = obj.q;
            obj.n_uvar = obj.p;
            
            obj.xprime = 0.5 * ones(1, obj.q);
            % bounds
            %init bound upper level
            obj.xu_bl = zeros(1, obj.p);
            obj.xu_bu = ones(1, obj.p);
           
           
            % init bound lower level
            obj.xl_bl = zeros(1, obj.q) ;
            obj.xl_bu = ones(1, obj.q) * 1;      
        end
        
        function [f, c] = evaluate_u(obj, xu, xl) 
            %-obj
            f = []
            c = [];
            
        end
        
        
        function [f, c] = evaluate_l(obj, xu, xl)
           
            %-obj
            f = 0;
            n = 1;
           for i = 1: obj. q
               index1  = xl(:, i) <  0.4696;
               index2  = xl(:, i) >= 0.4696 & xl(:, i) <= 0.5304;
               index3  = xl(:, 1) > 0.5304;
               
               fb1  =  -0.5 * exp(-0.5 * ((xl(:, i) - 0.4).^2 ./ (0.05^2)));
               fb2  =  -0.6 * exp(-0.5 * ((xl(:, i) - 0.5).^2 ./ (0.02^2)));
               fb3  =  -0.5 * exp(-0.5 * ((xl(:, i) - 0.6).^2 ./ (0.05^2)));
               
               fb1(index2) = 0; fb1(index3) = 0;
               fb2(index1) = 0; fb2(index3) = 0;
               fb3(index1) = 0; fb3(index2) = 0;
               f = f+ fb1 + fb2 + fb3;
           end
           
           f = f ./ n;
           
            %-con
            c = [];
         end
    end
end
