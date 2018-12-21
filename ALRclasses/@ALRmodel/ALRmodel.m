classdef ALRmodel < AutoLogistic 
    %ALRmodel Abstract class extending the AutoLogistic class to allow for regression.
    %   TODO: Detailed explanation goes here
    
    %TODO:
    % *
    
    properties (Dependent = true)
        X
        Beta
        Alpha
    end
    properties (SetAccess = private)
        P = 0;                                                %-Number of predictors.
    end
    properties (Access = private)
        PrivateX
        PrivateBeta
    end
    
    methods
        function obj = set.X(obj,V)
            %TODO: input checking
            obj.PrivateX = V;
            obj.N = size(V,1);
            obj.M = size(V,3);
            obj.P = size(V,2);
            obj.DimensionsOK = CheckDims(obj);
        end
        function X = get.X(obj)
            X = obj.PrivateX;
        end
        
        function obj = set.Beta(obj,V)
            %TODO: input checking
            if (isrow(V))
                obj.PrivateBeta = V';
            else
                obj.PrivateBeta = V;
            end
            obj.P = length(V);
            obj.DimensionsOK = CheckDims(obj);
        end
        function Beta = get.Beta(obj)
            Beta = obj.PrivateBeta;
        end
        
        function obj = set.Alpha(obj,~)
            warning('You cannot set Alpha; it is calculated from X and Beta')
        end
        function Alpha = get.Alpha(obj)
            [d1, d2, d3] = size(obj.X);
            if d2==length(obj.Beta)
                Alpha = zeros(d1,d3);
                for i=1:d3
                    Alpha(:,i) = obj.X(:,:,i)*obj.Beta;
                end
            else
                error('Can''t get Alpha: dimensions inconsistent.')
            end
        end
        
        function out = Predict(obj,newX)
            %TODO: jazz this up. 
            newobj = obj;
            newobj.X = newX;
            out = newobj.MarginalProbability();
        end
    end
    
    methods (Access=protected)
        function tf = CheckDims(obj)
            ShouldBeN = [size(obj.X,1), size(obj.Y,1), obj.Graph.numnodes];
            ShouldBeM = [size(obj.X,3), size(obj.Y,2)];
            ShouldBeP = [size(obj.X,2), length(obj.Beta)];
            tf = all(ShouldBeN==obj.N) && all(ShouldBeP==obj.P) ...
                && all(ShouldBeM==obj.M);
        end
    end
    
end

