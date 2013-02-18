require "./lib/Snippet"
require "bio"
require "amatch"

class ExperimentError < StandardError ; end
class RubyPrimerError < StandardError ; end


class Experiment

    attr_reader :mutatedTemplate, :PPsnippet, :DefaultPPtm, :forwardPrimer, 
    :forwardPrimerTemplate, :reversePrimer, :reversePrimerTemplate, :template,
    :experimentString, :ppTM

    @@AAcodes = 'acdefghiklmnpqrstvwy'
    @@defaultPPtm = 50.0

    def initialize(experimentString, template)
        @experimentString = experimentString
        @template = Snippet.new(template, template)
        @ppTM = @@defaultPPtm
        while not @forwardPrimer or @ppTM <= @@defaultPPtm or @forwardPrimer.length < 30
            @ppTM += 1
            self.setMutatedTemplate
            self.getForwardPrimer
            self.getReversePrimer
        end
        if ! self.mutatedTemplate.include?(@forwardPrimer.to_s) or ! self.mutatedTemplate.include?(@reversePrimer.reverse_complement)
            raise RubyPrimerError("Internal Error: Contact Simon")
        end
    end

    def setPPRegion(sequence, tm=nil)
        @PPsnippet = Snippet.new(sequence, @mutatedTemplate)
        if tm
            @PPsnippet = @PPsnippet.adjustTM(tm, ends=:both)
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
        @reversePrimer = Bio::Sequence::NA.new(@reversePrimerTemplate.snippet + @PPsnippet.snippet).reverse_complement
    end

    def printForwardPrimer
        "<span class='highlightPP'>#{@PPsnippet.to_s}</span>#{@forwardPrimerTemplate.to_s}"
    end

    def printReversePrimer
        rpPP = @PPsnippet.reverse_complement
        rpPT = @reversePrimerTemplate.reverse_complement
        "<span class='highlightPP'>#{rpPP}</span>#{rpPT}"
    end

    def log
        entry = "#{self.class}\n#{@template}, #{@experimentString}\nFP: #{@forwardPrimer.to_s}, RP: #{reversePrimer.to_s}\nPPtm: #{@PPsnippet.tm.round(1)}, FPT: #{@forwardPrimerTemplate.tm.round(1)}, RPT: #{@reversePrimerTemplate.tm.round(1)}\nPP Seq: #{@PPsnippet.to_s}"
        if @insertionSeq
            entry << "\nInsertion: #{@insertionSeq}"
        end
        return entry
    end

    def pcrConditions
        Hash[:cycles, 18, :annealing, [@forwardPrimerTemplate.tm, @reversePrimerTemplate.tm].min.round(1) - 5.0, :finalannealing, (@PPsnippet.tm - 5.0).round(1)]
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
        @PPsnippet = Snippet.new(fivePrimeSnippet.to_s + threePrimeSnippet.to_s, @mutatedTemplate)
        @PPsnippet = @PPsnippet.adjustTM(@ppTM, ends=:both)
        #self.setPPRegion(fivePrimeSnippet.to_s + threePrimeSnippet.to_s, tm=@ppTM)
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

    attr_reader :errorString

    def initialize(experimentString, template, errorString=nil)
        @errorString = errorString
        @experimentString = experimentString
        @template = template
    end

    def printData
        puts "Error Experiment Created"
        puts @experimentString
        puts @template
    end

    def log
        return "Error in Experiment (#{@errorString})\n#{@template}, #{experimentString}"
    end

end