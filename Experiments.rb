require "./Snippet"

class Experiment

	attr_reader :mutatedTemplate, :PPsnippet, :DefaultPPtm

	@@AAcodes = 'acdefghiklmnpqrstvwy'
	@@DefaultPPtm = 48.0

	def initialize(experimentString, template)
		@experimentString = experimentString
		@template = Snippet.new(template, template)
		@mutatedTemplate = self.getmutatedTemplate
		if @PPsnippet
			@PPsnippet = @PPsnippet.adjustTM(@@DefaultPPtm)
		end
	end
end



class DeletionExperiment < Experiment

	@@ExperimentRegex = /\A[#{@@AAcodes}]+-[#{@@AAcodes}]+\z/i

	def getmutatedTemplate
		templatebits = @experimentString.match /(?<fivePrime>[a-z]+)-(?<threePrime>[a-z]+)/i
		fivePrimeSnippet = @template.backTranslate(templatebits["fivePrime"])
		threePrimeSnippet = @template.backTranslate(templatebits["threePrime"])
		mutatedTemplate = @template.to_s[0..fivePrimeSnippet.end] + @template.to_s[threePrimeSnippet.start..-1]
		@PPsnippet = Snippet.new(fivePrimeSnippet.to_s + threePrimeSnippet.to_s, mutatedTemplate)
		return mutatedTemplate
	end

	def getForwardPrimer
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

@template = 'atgcccgctgaaacgaccgtatccggcgcgcaccccgccgccaaactgccgatttacatc'
@aaSeq = 'MPAETTVSGAHPAAKLPIYI'
@expString = 'MPAE-PIYI'
@expectedMutTemplate = 'atgcccgctgaaccgatttacatc'
@ExperimentObject = DeletionExperiment.new(@expString, @template)