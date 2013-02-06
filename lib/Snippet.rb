require 'Bio'

class SnippetError < StandardError ; end
class DNAFormatError < StandardError ; end
class DNAIndexError < StandardError ; end


class Bio::Sequence::NA

    def gc_content
        count = self.composition
        at = count['a'] + count['t'] + count['u']
        gc = count['g'] + count['c']
        0.0 if at + gc == 0
        return gc.quo(at + gc).to_f
    end

    def tm(ion=0.05, mismatch=0)
        mismatch_pc = mismatch.to_f / self.length * 100
        t = 81.5 + 
        16.6 * Math.log10(ion/(1 + 0.7 * ion)) + 
        0.41 * self.gc_percent - 
        500.0/self.length - mismatch_pc
        return t
    end
end


class Snippet

    attr_reader :snippet, :template, :start, :end, :BRSnippet, :BRTemplate

    def initialize(snippetSequence, templateSequence, start=nil, finish=nil)
        @template = templateSequence.strip.delete("\n\t")
        if start and finish
            @start = start
            @end = finish
            @snippet = @template[@start..@end]
        else
            @snippet = snippetSequence.strip.delete("\n\t")
            raise SnippetError.new("0 or >1 snippet occurrence found") if @template.scan(@snippet).length != 1
            @start = @template.index(@snippet)
            raise SnippetError unless @start 
            @end = @start + @snippet.length - 1
        end
        @BRSnippet = Bio::Sequence::NA.new(@snippet)
        @BRTemplate = Bio::Sequence::NA.new(@template)
        if @BRSnippet.illegal_bases.any? or @BRTemplate.illegal_bases.any? then raise DNAFormatError end        
    end        

    def method_missing(method_name, *args, &block)
        @BRSnippet.send(method_name, *args, &block)
    end

    def adjustTM(tm, ends=:both)
        newSnippet = self.clone
        if ends == :both
            while newSnippet.tm < tm
                newSnippet.start -= 1
                newSnippet.end += 1 if newSnippet.tm < tm
            end
            while newSnippet.tm > tm
                newSnippet.start += 1
                newSnippet.end -= 1 if newSnippet.tm > tm
            end
            return newSnippet
        end

        if ends == :left
            newSnippet.start +=1 while newSnippet.tm > tm
            newSnippet.start -=1 while newSnippet.tm < tm
            return newSnippet
        end

        newSnippet.end -=1 while newSnippet.tm > tm
        newSnippet.end +=1 while newSnippet.tm < tm
        return newSnippet
    end

    def backTranslate(substring)
        aaTrans = @BRTemplate.translate
        subStart = aaTrans.index(substring.upcase)
        raise SnippetError.new("Snippet sequence not found in template") unless subStart
        subStart = subStart * 3
        subEnd = subStart + substring.length * 3
        return Snippet.new(@template[subStart, substring.length * 3], @template)
    end

    def inspect
        {:sequence => @snippet, :start => @start, :end => @end, :id => self.object_id}
    end

    def to_s
        @snippet
    end

    def shift(numberOfBases)
        initialize(snippetSequence=nil, templateSequence=@template, start=@start+=numberOfBases, finish=@end+=numberOfBases)
    end

    def extend5(numberOfBases)
        self.start += (numberOfBases * -1)
    end

    def extend3(numberOfBases)
        self.end += numberOfBases
    end

    def start=(start)
        raise DNAIndexError.new("can't have a negative start!") if start < 0
        initialize(snippetSequence=nil, templateSequence=@template, start=start, finish=self.end)
    end

    def end=(endVal)
        raise DNAIndexError.new("end value > template length") if endVal >= @template.length
        initialize(snippetSequence=nil, templateSequence=@template, start=self.start, finish=endVal)
    end

    def length
        @snippet.length
    end

    def +(otherSnippet)
        snips = [self, otherSnippet].sort_by {|i| i.start}
        raise SnippetError.new("Snippets are not consecutive") unless snips[1].start - snips[0].end == 1
        return Snippet.new(snippetSequence = nil, templateSequence = @template, snippetSequence = snips[0].start, finish = snips[1].end)
    end
end
