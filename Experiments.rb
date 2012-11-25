require "./Snippet"

class ExperimentError < StandardError ; end

class Experiment

	attr_reader :mutatedTemplate, :PPsnippet, :DefaultPPtm, :forwardPrimer, 
	:forwardPrimerTemplate, :reversePrimer, :reversePrimerTemplate

	@@AAcodes = 'acdefghiklmnpqrstvwy'
	@@defaultPPtm = 50.0

	def initialize(experimentString, template)
		@experimentString = experimentString
		@template = Snippet.new(template, template)
		@@defaultPPtm -= 1
		while not @forwardPrimer or @forwardPrimer.length <= 30
			@@defaultPPtm += 1
			@mutatedTemplate = self.getmutatedTemplate
			self.getForwardPrimer
			@reversePrimer = self.getReversePrimer
		end
	end
end



class DeletionExperiment < Experiment

	@@ExperimentRegex = /\A[#{@@AAcodes}]+-[#{@@AAcodes}]+\z/i

	def getmutatedTemplate
		templatebits = @experimentString.match /(?<fivePrime>[a-z]+)-(?<threePrime>[a-z]+)/i
		fivePrimeSnippet = @template.backTranslate(templatebits["fivePrime"])
		threePrimeSnippet = @template.backTranslate(templatebits["threePrime"])
		raise ExperimentError.new("C-terminal Sequence comes before N-terminal Sequence") if fivePrimeSnippet.start > threePrimeSnippet.start
		mutatedTemplate = @template.to_s[0..fivePrimeSnippet.end] + @template.to_s[threePrimeSnippet.start..-1]
		@PPsnippet = Snippet.new(fivePrimeSnippet.to_s + threePrimeSnippet.to_s, mutatedTemplate)
		@PPsnippet = @PPsnippet.adjustTM(@@defaultPPtm, ends=:both)
		return mutatedTemplate
	end

	def getForwardPrimer
		@forwardPrimerTemplate = Snippet.new(snippetSequence = nil, templateSequence=@mutatedTemplate, start=(@PPsnippet.end + 1), finish=(@PPsnippet.end + 3))
		@forwardPrimerTemplate = @forwardPrimerTemplate.adjustTM(@PPsnippet.tm + 5.0, ends=:right)
		@forwardPrimer = Snippet.new(snippetSequence = nil, templateSequence=@mutatedTemplate, start=@PPsnippet.start, finish=@forwardPrimerTemplate.end)
	end
	
	def getReversePrimer
		@reversePrimerTemplate = Snippet.new(snippetSequence = nil, templateSequence=@mutatedTemplate, start = (@PPsnippet.start - 4), finish=(@PPsnippet.start - 1))
		@reversePrimerTemplate = @reversePrimerTemplate.adjustTM(@PPsnippet.tm + 5.0, ends=:left)
		reversePrimer = Bio::Sequence::NA.new(@reversePrimerTemplate.snippet + @PPsnippet.snippet)
		return reversePrimer.reverse_complement
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
