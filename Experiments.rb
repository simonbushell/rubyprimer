require "Snippet"

class Experiment
	def initialize(experimentString, template)
	end
end


class DeletionExperiment < Experiment

	@@ExperimentRegex = /\A[acdefghiklmnpqrstvwy]+-[acdefghiklmnpqrstvwy]+\z/i
end


class InsertionExperiment < Experiment
	
	@@ExperimentRegex = /\A[acdefghiklmnpqrstvwy]+[ACDEFGHIKLMNPRSTVWY]+[acdefghiklmnpqrstvwy]+\z/ 
end


class SubstituionExperiment < Experiment
	
	@@ExperimentRegex = /\A[acdefghiklmnpqrstvwy]+\z/ 
end