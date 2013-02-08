require "./lib/Snippet"
require "Bio"
require "amatch"

class ExperimentError < StandardError ; end


class Experiment

    attr_reader :mutatedTemplate, :PPsnippet, :DefaultPPtm, :forwardPrimer, 
    :forwardPrimerTemplate, :reversePrimer, :reversePrimerTemplate, :template,
    :experimentString

    @@AAcodes = 'acdefghiklmnpqrstvwy'
    @@defaultPPtm = 50.0

    def initialize(experimentString, template)
        @experimentString = experimentString
        @template = Snippet.new(template, template)
        @ppTM = @@defaultPPtm - 1
        while not @forwardPrimer or @forwardPrimer.length <= 30
            @ppTM += 1
            self.setMutatedTemplate
            self.getForwardPrimer
            @reversePrimer = self.getReversePrimer
        end
    end

    def setPPRegion(sequence, tm=nil)
        @PPsnippet = Snippet.new(sequence, @mutatedTemplate)
        if tm
            @PPsnippet = @PPsnippet.adjustTM(@ppTM, ends=:both)#
        end
    end

    def getBestCodon(codon, aminoacid)
        # include Amatch
        codonTable = Bio::CodonTable[1]
        possibleCodons = codonTable.revtrans(aminoacid)
        codonScores = Hash.new
        m = Amatch::Sellers.new(codon)
        possibleCodons.each {|x| codonScores[x] = m.match(x)}
        return codonScores.sort_by{|codon, score| score}[0][0]
    end

    def printData
        puts @experimentString
        puts "Forward Primer: #{@forwardPrimer.to_s} (#{@forwardPrimer.length} bp)"
        puts "Reverse Primer: #{@reversePrimer.to_s} (#{@reversePrimer.length} bp)"
        puts "\nPP tm: #{@PPsnippet.tm.round(1)}"
        puts "Forward PT tm: #{@forwardPrimerTemplate.tm.round(1)}"
        puts "Reverse PT tm: #{@reversePrimerTemplate.tm.round(1)}"
        puts "\nPP Sequence: #{@PPsnippet.to_s}"
        puts "rev PPSequence: #{@PPsnippet.reverse_complement}"
    end

    def getForwardPrimer
        # Calculates the required forward primer for the experiment and returns
        # it as a Snippet
        @forwardPrimerTemplate = Snippet.new(snippetSequence = nil, templateSequence=@mutatedTemplate, start=(@PPsnippet.end + 1), finish=(@PPsnippet.end + 3))
        @forwardPrimerTemplate = @forwardPrimerTemplate.adjustTM(@PPsnippet.tm + 5.0, ends=:right)
        @forwardPrimer = Snippet.new(snippetSequence = nil, templateSequence=@mutatedTemplate, start=@PPsnippet.start, finish=@forwardPrimerTemplate.end)
    end
    
    def getReversePrimer
        @reversePrimerTemplate = Snippet.new(snippetSequence = nil, templateSequence=@mutatedTemplate, start = (@PPsnippet.start - 4), finish=(@PPsnippet.start - 1))
        @reversePrimerTemplate = @reversePrimerTemplate.adjustTM(@PPsnippet.tm + 5.0, ends=:left)
        reversePrimer = Bio::Sequence::NA.new(@reversePrimerTemplate.snippet + @PPsnippet.snippet)
        reversePrimer.reverse_complement
    end
end


class DeletionExperiment < Experiment

    @@ExperimentRegex = /\A[#{@@AAcodes}]+-[#{@@AAcodes}]+\z/i

    def setMutatedTemplate
        # Returns the mutated template as a string, and calculates the required
        # primer-primer overlap region as a snippet of the mutated template

        templatebits = @experimentString.match /(?<fivePrime>[a-z]+)-(?<threePrime>[a-z]+)/i
        fivePrimeSnippet = @template.backTranslate(templatebits["fivePrime"])
        threePrimeSnippet = @template.backTranslate(templatebits["threePrime"])
        raise ExperimentError.new("C-terminal Sequence comes before N-terminal Sequence") if fivePrimeSnippet.start > threePrimeSnippet.start
        @mutatedTemplate = @template.to_s[0..fivePrimeSnippet.end] + @template.to_s[threePrimeSnippet.start..-1]
        self.setPPRegion(fivePrimeSnippet.to_s + threePrimeSnippet.to_s, tm=@ppTM)
    end
end


class InsertionExperiment < Experiment
    
    @@ExperimentRegex = /\A[#{@@AAcodes}]+[#{@@AAcodes}]+[#{@@AAcodes}]+\z/i
    @@defaultPPtm = 45.0

    attr_reader :insertionSeq

	def generateInsertion(proteinSequence)
	    optimisedCodons = Hash["G", "ggc", "E", "gaa", "D", "gat", "V", "gtg", 
	        "A", "gcg", "R", "cgc", "K", "aaa", "N", "aac", "M", "atg", "I",
	        "att", "T", "acc", "W", "tgg", "C", "tgc", "X", "taa", "Y", "tat",
	        "F", "ttt", "S", "agc", "Q", "cag", "H", "cat", "L", "ctg", "P", "ccg"]
	    generatedSequence = ""
	    proteinSequence.each_char { |x| generatedSequence << optimisedCodons[x]}
	    return generatedSequence
	end

	def getForwardPrimer
		insertSnippet = Snippet.new(@insertionSeq, @mutatedTemplate)
		@forwardPrimerTemplate = Snippet.new(snippetSequence = nil, templateSequence=@mutatedTemplate, start=(insertSnippet.end + 1), finish=(insertSnippet.end + 3))
		@forwardPrimerTemplate = @forwardPrimerTemplate.adjustTM(@PPsnippet.tm + 5.0, ends=:right)
		@forwardPrimer = Snippet.new(snippetSequence = nil, templateSequence=@mutatedTemplate, start=@PPsnippet.start, finish=@forwardPrimerTemplate.end)
	end

	def getReversePrimer
		super
		reversePrimer = Bio::Sequence::NA.new(@reversePrimerTemplate.snippet + @PPsnippet.snippet + @insertionSeq)
        reversePrimer.reverse_complement
	end

	def setMutatedTemplate
		preMutSnippet = @template.backTranslate(@experimentString.split('+')[0])
        @insertionSeq = self.generateInsertion(@experimentString.split('+')[1])
        postMutSnippet = @template.backTranslate(@experimentString.split('+')[2])

        @mutatedTemplate = @template.to_s[0..preMutSnippet.end] + @insertionSeq + @template.to_s[postMutSnippet.start .. -1]
        @PPsnippet = Snippet.new(preMutSnippet.to_s, @mutatedTemplate)
        @PPsnippet = @PPsnippet.adjustTM(@ppTM, ends=:left)
    end

    def printData
    	super
    	puts "Insertion Sequence: #{@insertionSeq}"
    end
end


class SubstitutionExperiment < Experiment
    
    @@ExperimentRegex = /\A(?<presub>[#{@@AAcodes}]+)[*](?<sub>[#{@@AAcodes}])[*](?<postsub>[#{@@AAcodes}])+\z/i

    def setMutatedTemplate
        preMutSnippet = @template.backTranslate(@experimentString.split('*')[0])
        mutationAA = @experimentString.split('*')[1]
        postMutSnippet = @template.backTranslate(@experimentString.split('*')[2])

        codonToChange = @template.to_s[preMutSnippet.end + 1 .. postMutSnippet.start - 1]
        raise ExperimentError.new("This requires more than a single codon change") unless codonToChange.length == 3

        @mutatedTemplate = @template.to_s[0..preMutSnippet.end] + self.getBestCodon(codonToChange, mutationAA) + @template.to_s[postMutSnippet.start .. -1]
        mutatedRegion = preMutSnippet.to_s + self.getBestCodon(codonToChange, mutationAA) + postMutSnippet.to_s
        self.setPPRegion(mutatedRegion, tm=@ppTM)
    end
end


class ErrorExperiment < Experiment

    def initialize(experimentString, template)
    end

end