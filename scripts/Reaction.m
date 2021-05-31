classdef Reaction < handle
    
    %REACTION creates a sealed class of an individual reaction, which contains information of reactants, stoichiometry and kinetic equations.
    
    properties        
        % A cell array of names of all reactants involved in this reaction.
        reactants
        %
        % State-change vector of this reaction.
        stateChangeVec
        %
        % Other factors that regulates the rate of this reaction but their
        % quantity does not change after the reaction occurs.
        regulators        
        %
        rateFunc
        %
        % e.g. For a reaction [2H2 + O2 -> 2H20]:
        %   reactantName: {'H2','O2','H2O'};
        %   stateChangeVec: [-2, -1, 2];
        % where the two properties match to each other.
    end
    
    methods        
        function obj = Reaction(reactantName, stateChangeVec, regulators, rateFunction)            
            %REACTION is the constructor function for class REACTION.
            % Input parameters for defining two properties should be
            % matched to each other.
            %
            if length(reactantName) ~= length(stateChangeVec)
                error('Number of reactants is not consistent with number of state-change vectors!')
            end
            
            obj.reactants = reactantName;
            obj.stateChangeVec = stateChangeVec; 
            obj.regulators = regulators;
            obj.rateFunc = rateFunction;
        end             
    end    
    
end