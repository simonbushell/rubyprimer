require "./Snippet"

class Experiment

	attr_reader :mutatedTemplate

	@@AAcodes = 'acdefghiklmnpqrstvwy'

	def initialize(experimentString, template)
		@experimentString = experimentString
		@template = Snippet.new(template, template)
		@mutatedTemplate = self.getmutateTemplate
	end
end



class DeletionExperiment < Experiment

	@@ExperimentRegex = /\A[#{@@AAcodes}]+-[#{@@AAcodes}]+\z/i

	def getmutateTemplate
		templatebits = @experimentString.match /(?<fivePrime>[a-z]+)-(?<threePrime>[a-z]+)/i
		fivePrimeSnippet = @template.backTranslate(templatebits["fivePrime"])
		threePrimeSnippet = @template.backTranslate(templatebits["threePrime"])
		mutatedTemplate = @template.to_s[0..fivePrimeSnippet.end] + @template.to_s[threePrimeSnippet.start..-1]
		return mutatedTemplate
	end

end




class InsertionExperiment < Experiment
	
	@@ExperimentRegex = /\A[#{@@AAcodes}]+[ACDEFGHIKLMNPRSTVWY]+[#{@@AAcodes}]+\z/ 
end


class SubstituionExperiment < Experiment
	
	@@ExperimentRegex = /\A(?<presub>[#{@@AAcodes}]+)[*](?<sub>[#{@@AAcodes}])[*](?<postsub>[#{@@AAcodes}])+\z/i

	def mutateTemplate
		substrings = @@ExperimentRegex.match(experimentString)
		templateRegex = /#{substrings["presub"]}[#{@@AAcodes}]#{substrings["postsub"]}/
	end
end