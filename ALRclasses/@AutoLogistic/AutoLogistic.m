classdef (Abstract) AutoLogistic
    %AUTOLOGISTIC An abstract top-level class for autologistic models.
    %   The purposes of this class are:
    %   1. To hold all properties and methods that are certain to be common to all
    %      subclasses.
    %   2. To hold static methods that are potentially useful for subclasses and
    %      for working with these objects.
    %   3. To specify abstract members that must be implemented in subclasses, but
    %      that might have implementations varying from subclass to subclass.

    %NOTE ON USING PRIVATE AND DEPENDENT PROPERTIES:
    % Desired behaviour: Any time we set Y or Graph, the following happens:
    %   1. N is updated.
    %   2. DimensionsOK is updated.
    % This is advantageous over just making N, and DimensionsOK dependent methods,
    % because it avoids having to re-calculate them each time.
    % 
    % Based on MATLAB class guidelines, this means Y and Graph have to be 
    % dependent properties (to avoid initialization order dependence on
    % saving/loading). But of course these properties are not to be calculated on 
    % access like regular dependent properties.
    % 
    % So the solution is as follows:
    %   * Make new properties PrivateAlpha, PrivateGraph. These will hold 
    %     the actual objects. They are not dependent, and they have Access=private.
    %   * Make the original properties Alpha, Graph dependent.
    %   * Write get methods for Alpha, Graph, that just return the value of their
    %     private counterparts.
    %   * Write set methods for Alpha, Graph, that set values for the private
    %    counterparts (also updating N and DimensionOK as desired).
    
    %NOTE ON USING FUNCTION HANDLES IN PROPERTIES:
    % We can let lambda be a function like Lambda = @(i,j,obj) some_fcn(i,j,obj.X).
    % Then class methods can use this function, e.g.:
    %   function out = some_method(obj,i,j)
    %       out = obj.Lambda(i,j,obj);
    %   end
    % See the CEProblem.m class def for examples of this.
    
    %TODO:
    % [] **IMPORTANT** think about whether certain methods (like Mu() conditional
    %    probs, prob table, etc.) should actually be dependent properties!!!
    %    This would allow us to, e.g., subset like CZO.Mu(1:10).
    % [] Turn this into a package and move functions like MakeGraph into the package.
    % [x] Maybe implement a "configuration" property that has true/false, and ensure
    %    that Y is always numeric and matching the coding.  Otherwise need to
    %    frequently convert Y into numeric for computations.... or make a "cast"
    %    method... (***done already, need to check speed issues***)
    % [x] Implement a "ConditionalProbability" method that gives the conditional
    %    probability of each variable being "high". (or optionally to give conditional
    %    probability that Yi = yi, used for computing PSlik)
    % [] Do checking to make sure the whole set of code works right when M>1
    %   (replicates). Also, call them something other than "replicates?"
    
    %============= PROPERTIES =======================================================
    properties (Abstract)
        Alpha                %-The unary parameters (equals X*Beta in the ALR model).
                             % An N-by-M matrix if there are M replicates.
        AssociationMatrix    %-The pairwise parameters, a symmetric N-by-N matrix.
    end
    
    properties
        Description
        States = {'low'; 'high'};
        Centered
    end
    
    properties (Dependent = true)
        Y                        %-Responses, N-by-M matrix (M>1 only if replicates).
        Graph                              %-Undirected graph object with N vertices.
        Coding                              %-Numeric coding of the binary variables.
        Configuration      %-Coding-independent representation of Y using the states.
    end
    
    properties (SetAccess = protected)
        N = 0;                                    %-Number of variables in the graph.
        M = 0;                                                %-Number of replicates.
        DimensionsOK                    %-A flag set by dimension consistency checks.
    end
    
    properties (Access = private)
        PrivateY                                                     %-Used to set Y.
        PrivateGraph = graph();                                  %-Used to set Graph.
        PrivateCoding = [-1 1];                                 %-Used to set Coding.
    end
    
    %============= METHODS ==========================================================
    methods (Abstract, Access = protected)
        CheckDims(obj)    %-A function for checking consistency of all objects' dims.
    end
    
    methods (Abstract)
        Predict(obj)      %-Predictions for new data (implemented in subclasses where
                          % data is found.)
    end
    
    methods
        %--- Externally defined methods ---
        out = Mu(obj,nodes,reps)    %-Calculates the vector of centering adjustments.
        out = Negpotential(obj,varargin)      %-Calculates the negpotential function.
        out = ProbabilityTable(obj,force)                         %-Tabluate the PMF.
        out = MarginalProbability(obj,nodes)    %-Compute or estimate marginal probs.
        out = Transform(obj,coding,centered)              %-Parameter transformation.
        out = GibbsSample(obj,varargin)                              %-Gibbs sampler.
        out = PerfectSample(obj)                                   %-Perfect sampler.
        out = PerfectSample2(obj)
        out = ConditionalProbability(obj,nodes,reps)    %-Compute the marginal probs.
        out = PseudoLikelihood(obj)      %-Compute the negative log pseudolikelihood.
        %----------------------------------
        
        function obj = set.Y(obj,V)
            validateattributes(V,{'double','logical','categorical'},{'2d'})
            newV = obj.CheckY(V);
            obj.PrivateY = newV;
            obj.N = size(newV,1);
            obj.M = size(newV,2);
            obj.DimensionsOK = CheckDims(obj);
        end
        function Y = get.Y(obj)
            Y = obj.PrivateY;
        end
        
        function config = get.Configuration(obj)
            config = categorical(obj.Y,obj.Coding,obj.States);
        end
        
        function obj = set.Graph(obj,G)
            %TODO: input checking
            obj.PrivateGraph = G;
            obj.N = G.numnodes;
            obj.DimensionsOK = CheckDims(obj);
        end
        function Graph = get.Graph(obj)
            Graph = obj.PrivateGraph;
        end
        
        function obj = set.Coding(obj,v)
            validateattributes(v,{'double'},{'vector','numel',2,'increasing'})
            if iscolumn(v) 
                v = v';
            end
            newY = zeros(size(obj.Y));
            newY(obj.Y==obj.Coding(1)) = v(1);
            newY(obj.Y==obj.Coding(2)) = v(2);
            obj.PrivateY = newY;
            obj.PrivateCoding = v;
        end
        function cod = get.Coding(obj)
            cod = obj.PrivateCoding;
        end
                
        function config = MAP(obj)
            %TODO: Function to find the MAP estimate.
            %Check if this is even possible for all casses (using max flow
            %way), implement only if it is. (external function?)
        end
        
    end
    
    methods (Access = protected)
        %--- Externally defined methods ---
        newY = CheckY(obj,Y)      %Ensures that Y is set properly and matches coding.
        %----------------------------------
    end
    
    methods (Static)
        %--- Externally defined methods ---
        G = MakeGraph(r,c,k)                              %-Produces a regular graph.
        %----------------------------------
    end

end

