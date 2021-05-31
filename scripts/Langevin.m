classdef Langevin < handle
    
    %LANGEVIN Creates a reaction system that can evolve with given kinetic rate functions using based on solution to Chemical Langevin Equations (CLE).
    %
    properties
        %
        % A cell array of all reactions. See class REACTION.
        reactions
        %
        % Substance pool that contains number of substances involved in all
        % reactions.
        currSubstances
    end
    
    methods
        function obj = Langevin(reactions, initSubstances)
            
            %LANGEVIN is the contructor function
            %
            % Check the correct data type            
            if ~isa(reactions,'cell')
                error('Error. \nThe first input must be a containers.Map, not a %s', class(reactions))
            end
            
            for i = 1:length(reactions)
            if ~isa(reactions{i},'Reaction')
                error('Error. \nElements in the first input must be a Reaction, not a %s', class(reactions{i})) 
            end
            end
            
            if ~isa(initSubstances,'containers.Map')
                error('Error. \nSecond input must be a containers.Map, not a %s', class(initSubstances))
            end
            
            obj.reactions = reactions;
            obj.currSubstances = initSubstances; %containers.Map
        end
        
        function [substanceUpdate, reactionUpdate] = evolution(obj, dt, whiteNoise)
            %EVOLUTION Evolves the reaction system with defined rate functions.
            %   INPUT:
            %       dt: time leap
            %       whiteNoise: true: SDE; false: ODE
            %   OUTPUT:
            %       substanceUpdate: updated value for each substance (containers.Map)
            %       reactionUpdate:  updated value for each reaction
            %
            if nargin < 3, whiteNoise = true; end
            
            substancesVal_before = cell2mat(obj.currSubstances.values);
            reactionUpdate = zeros(1,length(obj.reactions));
            
            for i = 1:length(obj.reactions)
                reaction = obj.reactions{i};
                currRegulators = zeros(1,length(reaction.regulators));
                
                for j = 1:length(reaction.regulators)
                    currRegulators(j) = obj.currSubstances(reaction.regulators{j});
                end
                
                reactionRate = reaction.rateFunc(currRegulators);
                reactionUpdate(i) = reactionRate * dt;
                
                for j = 1:length(reaction.reactants)
                    try
                        obj.currSubstances(reaction.reactants{j}) = ...
                        obj.currSubstances(reaction.reactants{j}) + ...
                        reactionRate * reaction.stateChangeVec(j) * dt + ...
                         sqrt(reactionRate * dt * (reactionRate > 0)) * reaction.stateChangeVec(j) * randn * whiteNoise;
                    catch
                        % If reaction rate is less than zero, error would
                        % occur when adding white noise.
                        disp(reaction.reactants)
                        disp(reactionRate)
                        error('Reaction rate is less than zero!')
                    end
                    
                    % Moderate inaccuracy would occur in LANGEVIN
                    % simulation, and therefore the number of substances
                    % can be slightly less than zero. An artificial
                    % correction is added here.
                    if obj.currSubstances(reaction.reactants{j}) < 0
                        obj.currSubstances(reaction.reactants{j}) = 0;
                    end
                    
                end
            end
            updateVal = cell2mat(obj.currSubstances.values) - substancesVal_before;
            substanceUpdate = containers.Map(obj.currSubstances.keys, updateVal);
        end
        
    end
end