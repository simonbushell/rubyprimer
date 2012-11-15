require 'Snippet'
require 'Experiments'

class RubyPrimerError < StandardError ; end


class RubyPrimer

	@@AAcodes = 'acdefghiklmnpqrstvwy'

	def initialize(experiment, template)
		@template = template
		@experiment = self.createExperiment(experiment)
	end

	def createExperiment(experimentString)
		insertionRegex = /\A[acdefghiklmnpqrstvwy]+[ACDEFGHIKLMNPRSTVWY]+[acdefghiklmnpqrstvwy]+\z/ 
		substittionRegex = /\A[acdefghiklmnpqrstvwy]+\z/

		return DeletionExperiment.new(experimentString, self.template) if !!(experimentString =~ DeletionExperiment.ExperimentRegex)
		return InsertionExperiment.new(experimentString, self.template) if !!(experimentString =~ InsertionExperiment.ExperimentRegex) 
		return SubstitutionExperiment.new(experimentString, self.template) if !!(experimentString =~ SubstitutionExperiment.ExperimentRegex)

		raise RubyPrimerError.new("Experiment String in incorrect format")
	end
end