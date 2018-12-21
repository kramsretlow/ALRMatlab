classdef ALRadaptive < ALRmodel & matlab.mixin.CustomDisplay
    %ALRadaptive A class for autologistic regression models with adaptive smoothing.
    %   TODO: Detailed explanation goes here
    
    %TODO:
    % *
    
    properties
        Lambda              %-Lambda should be a handle to a function of (obj, i, j).
        Gamma                   %-Parameter vector to be used for adaptive smoothing.
    end
    
    properties (Dependent = true)
        AssociationMatrix
    end
    
    methods
        function obj = ALRadaptive(varargin)
            %Constructor method. Accepts parameter/value pairs for setting the public
            %methods of the class. NOTE: inputParser by default is case insensitive 
            %and supports partial matching. 
            
            %Parse the inputs.
            p = inputParser;
            p.FunctionName = 'ALRmodel constructor';
            addParameter(p,'Description','An autologistic regression model',@ischar)
            addParameter(p,'X',[],@isnumeric)
            addParameter(p,'Beta',[],@isvector)
            addParameter(p,'Gamma',[],@isvector)
            addParameter(p,'Coding',[-1 1], ...
                @(x) validateattributes(x,{'numeric'},{'size',[1 2],'increasing'}))
            addParameter(p,'Centered',false, ...
                @(x) validateattributes(x,{'logical'},{'size',[1 1]}))
            addParameter(p,'Y',[], ...
                @(x) validateattributes(x,{'double','logical','categorical'}, {'2d'}))
            addParameter(p,'Lambda',[], ...
                @(x) validateattributes(x,{'function_handle'},{'scalar'}))
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
            obj.Gamma = p.Results.Gamma;
        end
        
        function obj = set.Lambda(obj,f)
            %TODO: input checking. V should be a function handle.
            obj.Lambda = f;
        end
        
        function AM = get.AssociationMatrix(obj)
            error(['Association matrix for ALRadaptive hasn''t ' ...
                   'been implemented yet.'])
        end
    end
    
    methods (Access=protected)
        function propgrp = getPropertyGroups(~)
            %This code is stolen directly from MATLAB OOP documentation, for changing
            %the property order in disp.
            proplist = {'Description','Coding','Centered','Y','X','Beta','Gamma', ...
                        'Lambda', 'N', 'P', 'M','Graph','DimensionsOK'};
            propgrp = matlab.mixin.util.PropertyGroup(proplist);
        end
    end
end


