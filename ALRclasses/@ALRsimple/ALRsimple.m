classdef ALRsimple < ALRmodel & matlab.mixin.CustomDisplay
    %ALRmodel A class for autologistic regression models with simple smoothing.
    %   TODO: Detailed explanation goes here
    
    %TODO:
    % *
    
    properties
        Lambda = 0;
    end
    
    properties (Dependent = true)
        AssociationMatrix
    end
    
    methods
        function obj = ALRsimple(varargin)
            %Constructor method. Accepts parameter/value pairs for setting the public
            %methods of the class. NOTE: inputParser by default is case insensitive 
            %and supports partial matching. 
            
            %Parse the inputs.
            p = inputParser;
            p.FunctionName = 'ALRmodel constructor';
            addParameter(p,'Description','An autologistic regression model',@ischar)
            addParameter(p,'X',[],@isnumeric)
            addParameter(p,'Beta',[],@isvector)
            addParameter(p,'Coding',[-1 1], ...
                @(x) validateattributes(x,{'numeric'},{'size',[1 2],'increasing'}))
            addParameter(p,'Centered',false, ...
                @(x) validateattributes(x,{'logical'},{'size',[1 1]}))
            addParameter(p,'Y',[], ...
                @(x) validateattributes(x,{'double','logical','categorical'}, {'2d'}))
            addParameter(p,'Lambda',0, ...
                @(x) validateattributes(x,{'numeric'},{'scalar'}))
            addParameter(p,'Graph',graph(),@(x) isa(x,'graph'))
            parse(p,varargin{:})
            
            %Create the object.
            obj.Description = p.Results.Description;
            obj.Coding = p.Results.Coding;
            obj.Centered = p.Results.Centered;
            obj.Y = p.Results.Y;
            obj.Lambda = p.Results.Lambda;
            obj.Graph = p.Results.Graph;
            obj.X = p.Results.X;
            obj.Beta = p.Results.Beta;
        end
        
        function obj = set.Lambda(obj,V)
            %TODO: input checking. V should be a scalar.
            obj.Lambda = V;
        end
        
        function AM = get.AssociationMatrix(obj)
                AM = obj.Lambda*adjacency(obj.Graph);
        end
        
    end
    
    methods (Access=protected)
        function propgrp = getPropertyGroups(~)
            %This code is stolen directly from MATLAB OOP documentation, for changing
            %the property order in disp.
            proplist = {'Description','Coding','Centered','Y','X','Beta', ...
                        'Lambda', 'N', 'P', 'M','Graph','DimensionsOK'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
end





