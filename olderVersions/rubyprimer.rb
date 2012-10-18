require 'Bio'

class IncorrectParameterError < StandardError ; end

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



class PrimerDesignExperiment

	INITIAL_PP = 50.0

	attr_reader :templateDNA, :experiment, :templateProtein

	def initialize(template, experimentString=nil)
		@templateDNA = Bio::Sequence::NA.new(template)
		@experiment = experimentString
		@templateProtein = Bio::Sequence::AA.new(@templateDNA.translate)
	end

	def AAtoNA(substring, template=nil)
		template = self.templateDNA unless template
		subStart = self.templateProtein.index(substring)
		return nil unless subStart
		subStart = subStart * 3
		subEnd = subStart + substring.length * 3
		{:sequence => self.templateDNA[subStart, substring.length * 3], :start => subStart, :end => subEnd}
	end

	def adjustRegiontoTm(adjustedtm, seqHash, direction, tolerance=3.0)
		
	end
end


class DeletionExperiment < PrimerDesignExperiment

	attr_reader :ppRegion

	def findFlanks
		flank1, flank2 = self.experiment.upcase.split("-").collect {|x| self.AAtoNA(x)}
		@ppRegion = Bio::Sequence::NA.new(flank1[:sequence] + flank2[:sequence])
	end
end


template = "atgcccgctgaaacgaccgtatccggcgcgcaccccgccgccaaactgccgatttacatc
ctgccctgcttcctttggataggcatcgtcccctttaccttcgcgctcaaactgaaaccg
tcgcccgacttttaccacgatgccgccgccgcagccggcctgattgtcctgttgttcctc
acggcaggaaaaaaactgtttgatgtcaaaatccccgccatcagcttccttctgtttgca
atggcggcgttttggtatcttcaggcacgcctgatgaacctgatttaccccggtatgaac
gacatcgtctcttggattttcatcttgctcgccgtcagcgcgtgggcctgccggagcttg
gtcgcacacttcggacaagaacgcatcgtgaccctgtttgcctggtcgctgcttatcggc
tccctgcttcaatcctgcatcgtcgtcatccagtttgccggctgggaagacacccctctg
tttcaaaacatcatcgtttacagcgggcaaggcgtaatcggacacatcgggcagcgcaac
aacctcggacactacctcatgtggggcatactcgccgccgcctacctcaacggacaacga
aaaatccccgccgccctcggcgtaatctgcctgattatgcagaccgccgttttaggtttg
gtcaactcgcgcaccatcttgacctacatagccgccatcgccctcatccttcccttctgg
tatttccgttcggacaaatccaacaggcggacgatgctcggcatagccgcagccgtattc
cttaccgcgctgttccaattttccatgaacaccattctggaaacctttactggcatccgc
tacgaaactgccgtcgaacgcgtcgccaacggcggtttcacagacttgccgcgccaaatc
gaatggaataaagcccttgccgccttccagtccgccccgatattcgggcacggctggaac
agttttgcccaacaaaccttcctcatcaatgccgaacagcacaacatatacgacaacctc
ctcagcaacttgttcacccattcccacaacatcgtcctccaactccttgcagagatggga
atcagcggcacgcttctggttgccgcaaccctgctgacgggcattgccgggctgcttaaa
cgccccctgacccccgcatcgcttttcctaatctgcacgcttgccgtcagtatgtgccac
agtatgctcgaatatcctttgtggtatgtctatttcctcatccctttcggactgatgctc
ttcctgtcccccgcagaggcttcagacggcatcgccttcaaaaaagccgccaatctcggc
atactgaccgcctccgccgccatattcgcaggattgctgcacttggactggacatacacc
cggctggttaacgccttttcccccgccactgacgacagtgccaaaaccctcaaccggaaa
atcaacgagttgcgctatatttccgcaaacagtccgatgctgtccttttatgccgacttc
tccctcgtaaacttcgccctgccggaataccccgaaacccagacttgggcggaagaagca
accctcaaatcactaaaataccgcccccactccgccacctaccgcatcgccctctacctg
atgcggcaaggcaaagttgcagaagcaaaacaatggatgcgggcgacacagtcctattac
ccctacctgatgccccgatacgccgacgaaatccgcaaactgcccgtatgggcgccgctg
ctacccgaactgctcaaagactgcaaagccttcgccgccgcgcccggtcatccggaagca
aaaccctgcaaatga"

expString = "YFRSDKS-LTPAS"

exp = DeletionExperiment.new(template, expString)

prot = "MPAETTVSGAHPAAKLPIYILPCFLWIGIVPFTFALKLKPSPDFYHDAAAAA
GLIVLLFLTAGKKLFDVKIPAISFLLFAMAAFWYLQARLMNLIYPGMNDIVSWIFILLAV
SAWACRSLVAHFGQERIVTLFAWSLLIGSLLQSCIVVIQFAGWEDTPLFQNIIVYSGQGV
IGHIGQRNNLGHYLMWGILAAAYLNGQRKIPAALGVICLIMQTAVLGLVNSRTILTYIAA
IALILPFWYFRSDKSNRRTMLGIAAAVFLTALFQFSMNTILETFTGIRYETAVERVANGG
FTDLPRQIEWNKALAAFQSAPIFGHGWNSFAQQTFLINAEQHNIYDNLLSNLFTHSHNIV
LQLLAEMGISGTLLVAATLLTGIAGLLKRPLTPASLFLICTLAVSMCHSMLEYPLWYVYF
LIPFGLMLFLSPAEASDGIAFKKAANLGILTASAAIFAGLLHLDWTYTRLVNAFSPATDD
SAKTLNRKINELRYISANSPMLSFYADFSLVNFALPEYPETQTWAEEATLKSLKYRPHSA
TYRIALYLMRQGKVAEAKQWMRATQSYYPYLMPRYADEIRKLPVWAPLLPELLKDCKAFA
AAPGHPEAKPCK*"
# puts exp.templateProtein

# puts exp.AAtoNA('DKSNRR').inspect

# puts exp.findFlanks

#puts exp.ppRegion

puts exp.templateProtein


 