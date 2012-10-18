require 'Bio'

class SnippetError < StandardError ; end
class DNAFormatError < StandardError ; end


class Bio::Sequence::NA

	def gc_content
	  count = self.composition
	  at = count['a'] + count['t'] + count['u']
	  gc = count['g'] + count['c']
	  return 0.0 if at + gc == 0
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

	def initialize(snippetSequence, templateSequence)
		@snippet = snippetSequence
		@template = templateSequence
		@BRSnippet = Bio::Sequence::NA.new(@snippet)
		@BRTemplate = Bio::Sequence::NA.new(@snippet)
		if @BRSnippet.illegal_bases.any? or @BRTemplate.illegal_bases.any? then raise DNAFormatError end
		@start = @template.index(@snippet)
		raise SnippetError unless @start 
		@end = @start + @snippet.length - 1 
	end		

	def method_missing(method_name, *args, &block)
		@BRSnippet.send(method_name, *args, &block)
	end

	def backTranslate(aaSeq)
	end

	def to_s
		{:sequence => @snippet, :start => @start, :end => @end}
	end

	def start=(start)
		initialize(@template[start..@end], @template)
	end

	def end=(endVal)
		initialize(@template[@start..endVal], @template)
	end
end


# @testtemplate = "atgcccgctgaaacgaccgtatccggcgcgcaccccgccgccaaactgccgatttacatc
# 		ctgccctgcttcctttggataggcatcgtcccctttaccttcgcgctcaaactgaaaccg
# 		tcgcccgacttttaccacgatgccgccgccgcagccggcctgattgtcctgttgttcctc
# 		acggcaggaaaaaaactgtttgatgtcaaaatccccgccatcagcttccttctgtttgca
# 		atggcggcgttttggtatcttcaggcacgcctgatgaacctgatttaccccggtatgaac
# 		gacatcgtctcttggattttcatcttgctcgccgtcagcgcgtgggcctgccggagcttg
# 		gtcgcacacttcggacaagaacgcatcgtgaccctgtttgcctggtcgctgcttatcggc
# 		tccctgcttcaatcctgcatcgtcgtcatccagtttgccggctgggaagacacccctctg
# 		tttcaaaacatcatcgtttacagcgggcaaggcgtaatcggacacatcgggcagcgcaac
# 		aacctcggacactacctcatgtggggcatactcgccgccgcctacctcaacggacaacga
# 		aaaatccccgccgccctcggcgtaatctgcctgattatgcagaccgccgttttaggtttg
# 		gtcaactcgcgcaccatcttgacctacatagccgccatcgccctcatccttcccttctgg
# 		tatttccgttcggacaaatccaacaggcggacgatgctcggcatagccgcagccgtattc
# 		cttaccgcgctgttccaattttccatgaacaccattctggaaacctttactggcatccgc
# 		tacgaaactgccgtcgaacgcgtcgccaacggcggtttcacagacttgccgcgccaaatc
# 		gaatggaataaagcccttgccgccttccagtccgccccgatattcgggcacggctggaac
# 		agttttgcccaacaaaccttcctcatcaatgccgaacagcacaacatatacgacaacctc
# 		ctcagcaacttgttcacccattcccacaacatcgtcctccaactccttgcagagatggga
# 		atcagcggcacgcttctggttgccgcaaccctgctgacgggcattgccgggctgcttaaa
# 		cgccccctgacccccgcatcgcttttcctaatctgcacgcttgccgtcagtatgtgccac
# 		agtatgctcgaatatcctttgtggtatgtctatttcctcatccctttcggactgatgctc
# 		ttcctgtcccccgcagaggcttcagacggcatcgccttcaaaaaagccgccaatctcggc
# 		atactgaccgcctccgccgccatattcgcaggattgctgcacttggactggacatacacc
# 		cggctggttaacgccttttcccccgccactgacgacagtgccaaaaccctcaaccggaaa
# 		atcaacgagttgcgctatatttccgcaaacagtccgatgctgtccttttatgccgacttc
# 		tccctcgtaaacttcgccctgccggaataccccgaaacccagacttgggcggaagaagca
# 		accctcaaatcactaaaataccgcccccactccgccacctaccgcatcgccctctacctg
# 		atgcggcaaggcaaagttgcagaagcaaaacaatggatgcgggcgacacagtcctattac
# 		ccctacctgatgccccgatacgccgacgaaatccgcaaactgcccgtatgggcgccgctg
# 		ctacccgaactgctcaaagactgcaaagccttcgccgccgcgcccggtcatccggaagca
# 		aaaccctgcaaatga"

# 		@snippetseq = "gtcgaacgcgtcgccaacggcggtttcacaga"

# 		@snippet = Snippet.new(@snippetseq, @testtemplate)
