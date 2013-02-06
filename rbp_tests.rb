

require "./lib/Snippet"
require "./lib/Experiments"
require "test/unit"
require "Bio"


class TestSnippet < Test::Unit::TestCase 

	def setup
		@testtemplate = "atgcccgctgaaacgaccgtatccggcgcgcaccccgccgccaaactgccgatttacatc
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

		@snippetseq = "gtcgaacgcgtcgccaacggcggtttcacaga"

		@snippet = Snippet.new(@snippetseq, @testtemplate)
	end

	def test_initialize
		assert(@snippet)
	end

	def test_snippet_validation
		assert_raise(SnippetError) {snip2 = Snippet.new('agagagagagagagagagaga', @testtemplate)}
		assert_raise(DNAFormatError) {snip2 = Snippet.new('gtcgar', 'tagtaggtcgartagtag')}
	end

	def test_tm
		assert_equal(@snippet.tm, 69.44989226541517)
	end

	def test_methods_missing
		assert(@snippet.gc_content)
		assert_send [@snippet, :gc_content]
		assert_send [@snippet, :illegal_bases]
	end

	def test_change_start
		template = 'tgatcgatcgattcgatcgatcgatcggctagct'
		snip = 'tcgatcgatc'
		snippet = Snippet.new(snip, template)
		init_snip = snippet.clone
		snippet.start = 0
		assert_equal(snippet.snippet, 'tgatcgatcgattcgatcgatc')
		assert_equal(snippet.start, 0)
		assert_equal(snippet.end, init_snip.end)
		assert_equal(snippet.template, init_snip.template)
		assert_not_equal(snippet.start, init_snip.start)
		assert_not_equal(snippet.BRSnippet, init_snip.BRSnippet)
	end

	def test_change_end
		template = 'tgatcgatcgattcgatcgatcgatcggctagct'
		snip = 'tcgatcgatc'
		snippet = Snippet.new(snip, template)
		init_snip = snippet.clone
		snippet.end = 26
		assert_equal(snippet.snippet, 'tcgatcgatcgatcg')
		assert_equal(snippet.end, 26)
		assert_equal(snippet.start, init_snip.start)
		assert_equal(snippet.template, init_snip.template)
		assert_not_equal(snippet.end, init_snip.end)
		assert_not_equal(snippet.BRSnippet, init_snip.BRSnippet)
	end

	def test_back_translate
		aaseq = "DKSNRR"
		d = @snippet.backTranslate(aaseq)
		assert_instance_of(Snippet, d)
		assert_equal(aaseq, d.translate)
		assert_equal(d.start, 732)
		assert_equal(d.end, 749)
	end

	def test_shift
		template = 'tgatcgatcgattcgatcgatcgatcggctagct'
		snip = 'tcgatcgatc'
		snippet = Snippet.new(snip, template)
		init_snip = snippet.clone
		snippet.shift(3)
		assert_equal(snippet.start - init_snip.start, 3)
		assert_equal(snippet.end - init_snip.end, 3)
		assert_equal(snippet.to_s, 'atcgatcgat')
		assert_equal(snippet.length, init_snip.length)
	end

	def test_adjustTM
		assert_equal(@snippet.tm, 69.44989226541517)
		d = @snippet.clone
		f = d.adjustTM(45)
		assert_not_equal(f.tm, 69.44989226541517)
		assert_not_equal(@snippet.snippet, f.snippet)
		assert(f.tm < 46)
	end

	def test_concatenation
		snip1 = Snippet.new('gacatcgtctctt', @testtemplate)
		snip2 = Snippet.new('ggattttcatctt', @testtemplate)
		snip2b = Snippet.new('gccttcaaaaaa', @testtemplate)
		snip3 = snip1 + snip2
		assert_equal(snip1.snippet + snip2.snippet, snip3.snippet)
		assert_raise(SnippetError) {snip4 = snip1 + snip2b}
		assert_equal(snip3.start, snip1.start)
		assert_equal(snip3.end, snip2.end)
	end

	def test_extend5
		test = @snippet.clone
		test.extend5(3)
		assert_equal(test.to_s, 'gcc'+@snippetseq)
		assert_equal(test.start, @snippet.start - 3)
		test = @snippet.clone
		test.extend5(-3)
		assert_equal(test.to_s, @snippetseq[3..-1])
		assert_equal(test.start, @snippet.start + 3)
	end

	def test_extend3
		test = @snippet.clone
		test.extend3(3)
		assert_equal(test.to_s, @snippetseq + 'ctt')
		assert_equal(test.end, @snippet.end + 3)
		test = @snippet.clone
		test.extend3(-3)
		assert_equal(test.to_s, @snippetseq[0..-4])
		assert_equal(test.end, @snippet.end - 3)
	end
end


class Test_DeletionExperiments < Test::Unit::TestCase 
	
	def setup
		@template = "atgcccgctgaaacgaccgtatccggcgcgcaccccgccgccaaactgccgatttacatc
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

		@aaSeq = 'MPAETTVSGAHPAAKLPIYILPCFLWIGIVPFTFALKLKPSPDFYHDAAAAAGLIVLLFL
		TAGKKLFDVKIPAISFLLFAMAAFWYLQARLMNLIYPGMNDIVSWIFILLAVSAWACRSL
		VAHFGQERIVTLFAWSLLIGSLLQSCIVVIQFAGWEDTPLFQNIIVYSGQGVIGHIGQRN
		NLGHYLMWGILAAAYLNGQRKIPAALGVICLIMQTAVLGLVNSRTILTYIAAIALILPFW
		YFRSDKSNRRTMLGIAAAVFLTALFQFSMNTILETFTGIRYETAVERVANGGFTDLPRQI
		EWNKALAAFQSAPIFGHGWNSFAQQTFLINAEQHNIYDNLLSNLFTHSHNIVLQLLAEMG
		ISGTLLVAATLLTGIAGLLKRPLTPASLFLICTLAVSMCHSMLEYPLWYVYFLIPFGLML
		FLSPAEASDGIAFKKAANLGILTASAAIFAGLLHLDWTYTRLVNAFSPATDDSAKTLNRK
		INELRYISANSPMLSFYADFSLVNFALPEYPETQTWAEEATLKSLKYRPHSATYRIALYL
		MRQGKVAEAKQWMRATQSYYPYLMPRYADEIRKLPVWAPLLPELLKDCKAFAAAPGHPEA
		KPCK-'

		@expString = 'DKSNRR-KRPLTP'
		@expectedMutTemplate = "atgcccgctgaaacgaccgtatccggcgcgcaccccgccgccaaactgccgatttacatc
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
		tatttccgttcggacaaatccaacaggcggaaacgccccctgacccccgcatcgcttttc
		ctaatctgcacgcttgccgtcagtatgtgccacagtatgctcgaatatcctttgtggtat
		gtctatttcctcatccctttcggactgatgctcttcctgtcccccgcagaggcttcagac
		ggcatcgccttcaaaaaagccgccaatctcggcatactgaccgcctccgccgccatattc
		gcaggattgctgcacttggactggacatacacccggctggttaacgccttttcccccgcc
		actgacgacagtgccaaaaccctcaaccggaaaatcaacgagttgcgctatatttccgca
		aacagtccgatgctgtccttttatgccgacttctccctcgtaaacttcgccctgccggaa
		taccccgaaacccagacttgggcggaagaagcaaccctcaaatcactaaaataccgcccc
		cactccgccacctaccgcatcgccctctacctgatgcggcaaggcaaagttgcagaagca
		aaacaatggatgcgggcgacacagtcctattacccctacctgatgccccgatacgccgac
		gaaatccgcaaactgcccgtatgggcgccgctgctacccgaactgctcaaagactgcaaa
		gccttcgccgccgcgcccggtcatccggaagcaaaaccctgcaaatga"
		@ExperimentObject = DeletionExperiment.new(@expString, @template)
	end

	def test_createMutatedTemplate
		assert_equal(@ExperimentObject.mutatedTemplate, @expectedMutTemplate.strip.delete("\n\t"))
	end

	def test_PPtm
		assert @ExperimentObject.PPsnippet
		#assert_operator @ExperimentObject.PPsnippet.tm, :<, 50.0
	end

	def test_makeForwardPrimer
		assert @ExperimentObject.forwardPrimer
		  # p @ExperimentObject.PPsnippet.to_s
		assert (@ExperimentObject.forwardPrimerTemplate.tm - @ExperimentObject.PPsnippet.tm) > 5.0
		 # p @ExperimentObject.forwardPrimer
		 # @ExperimentObject.forwardPrimerTemplate.to_s
		   # p @ExperimentObject.forwardPrimerTemplate.tm
		   # p @ExperimentObject.PPsnippet.tm
		assert(@expectedMutTemplate.include?(@ExperimentObject.forwardPrimerTemplate.snippet))
	end

	def test_makeReversePrimer
		assert @ExperimentObject.reversePrimer
		# p @ExperimentObject.forwardPrimer.to_s
		 # p @ExperimentObject.reversePrimer
		# p @ExperimentObject.PPsnippet.snippet
		rpObject = Bio::Sequence::NA.new(@ExperimentObject.reversePrimer)
		assert(@ExperimentObject.mutatedTemplate.include?(rpObject.reverse_complement))
		#p rpObject.reverse_complement
		 # p @ExperimentObject.reversePrimerTemplate.tm
		assert(@ExperimentObject.reversePrimer.start_with?(@ExperimentObject.PPsnippet.reverse_complement))
	end

	def test_printOutput
		#p @ExperimentObject.printData
	end
end


class Test_SubstitutionExperiments < Test::Unit::TestCase

	def setup
		@template = "atgcccgctgaaacgaccgtatccggcgcgcaccccgccgccaaactgccgatttacatc
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

		@aaSeq = 'MPAETTVSGAHPAAKLPIYILPCFLWIGIVPFTFALKLKPSPDFYHDAAAAAGLIVLLFL
		TAGKKLFDVKIPAISFLLFAMAAFWYLQARLMNLIYPGMNDIVSWIFILLAVSAWACRSL
		VAHFGQERIVTLFAWSLLIGSLLQSCIVVIQFAGWEDTPLFQNIIVYSGQGVIGHIGQRN
		NLGHYLMWGILAAAYLNGQRKIPAALGVICLIMQTAVLGLVNSRTILTYIAAIALILPFW
		YFRSDKSNRRTMLGIAAAVFLTALFQFSMNTILETFTGIRYETAVERVANGGFTDLPRQI
		EWNKALAAFQSAPIFGHGWNSFAQQTFLINAEQHNIYDNLLSNLFTHSHNIVLQLLAEMG
		ISGTLLVAATLLTGIAGLLKRPLTPASLFLICTLAVSMCHSMLEYPLWYVYFLIPFGLML
		FLSPAEASDGIAFKKAANLGILTASAAIFAGLLHLDWTYTRLVNAFSPATDDSAKTLNRK
		INELRYISANSPMLSFYADFSLVNFALPEYPETQTWAEEATLKSLKYRPHSATYRIALYL
		MRQGKVAEAKQWMRATQSYYPYLMPRYADEIRKLPVWAPLLPELLKDCKAFAAAPGHPEA
		KPCK-'

		@expString = 'LGVI*A*LIMQ'
		@expectedMutTemplate = "atgcccgctgaaacgaccgtatccggcgcgcaccccgccgccaaactgccgatttacatc
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
		tatttccgttcggacaaatccaacaggcggaaacgccccctgacccccgcatcgcttttc
		ctaatctgcacgcttgccgtcagtatgtgccacagtatgctcgaatatcctttgtggtat
		gtctatttcctcatccctttcggactgatgctcttcctgtcccccgcagaggcttcagac
		ggcatcgccttcaaaaaagccgccaatctcggcatactgaccgcctccgccgccatattc
		gcaggattgctgcacttggactggacatacacccggctggttaacgccttttcccccgcc
		actgacgacagtgccaaaaccctcaaccggaaaatcaacgagttgcgctatatttccgca
		aacagtccgatgctgtccttttatgccgacttctccctcgtaaacttcgccctgccggaa
		taccccgaaacccagacttgggcggaagaagcaaccctcaaatcactaaaataccgcccc
		cactccgccacctaccgcatcgccctctacctgatgcggcaaggcaaagttgcagaagca
		aaacaatggatgcgggcgacacagtcctattacccctacctgatgccccgatacgccgac
		gaaatccgcaaactgcccgtatgggcgccgctgctacccgaactgctcaaagactgcaaa
		gccttcgccgccgcgcccggtcatccggaagcaaaaccctgcaaatga"
		@ExperimentObject = SubstitutionExperiment.new(@expString, @template)
	end

	def test_mutatedTemplate
		assert_not_equal(@template, @ExperimentObject.mutatedTemplate)
	end

	def test_primers
		assert(@ExperimentObject.mutatedTemplate.include?(@ExperimentObject.forwardPrimer.to_s))
		revPrimer = Bio::Sequence::NA.new(@ExperimentObject.reversePrimer)
		assert(@ExperimentObject.mutatedTemplate.include?(revPrimer.reverse_complement))
	end

	def test_errors
		errorString = "PAISF*A*FAM"
		assert_raise ExperimentError do
			exp = SubstitutionExperiment.new(errorString, @template)
		end
		errorString = "*D*"
		assert_raise SnippetError do
			exp = SubstitutionExperiment.new(errorString, @template)
		end
		errorString = "M*A*AETTVSG"
		assert_raise SnippetError do
			exp = SubstitutionExperiment.new(errorString, @template)
		end
	end
end

class Test_InsertionExperiments < Test::Unit::TestCase

	def setup
		@template = "atgcccgctgaaacgaccgtatccggcgcgcaccccgccgccaaactgccgatttacatc
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
		aaaccctgcaaatga".strip.delete("\n\t")

		@aaSeq = 'MPAETTVSGAHPAAKLPIYILPCFLWIGIVPFTFALKLKPSPDFYHDAAAAAGLIVLLFL
		TAGKKLFDVKIPAISFLLFAMAAFWYLQARLMNLIYPGMNDIVSWIFILLAVSAWACRSL
		VAHFGQERIVTLFAWSLLIGSLLQSCIVVIQFAGWEDTPLFQNIIVYSGQGVIGHIGQRN
		NLGHYLMWGILAAAYLNGQRKIPAALGVICLIMQTAVLGLVNSRTILTYIAAIALILPFW
		YFRSDKSNRRTMLGIAAAVFLTALFQFSMNTILETFTGIRYETAVERVANGGFTDLPRQI
		EWNKALAAFQSAPIFGHGWNSFAQQTFLINAEQHNIYDNLLSNLFTHSHNIVLQLLAEMG
		ISGTLLVAATLLTGIAGLLKRPLTPASLFLICTLAVSMCHSMLEYPLWYVYFLIPFGLML
		FLSPAEASDGIAFKKAANLGILTASAAIFAGLLHLDWTYTRLVNAFSPATDDSAKTLNRK
		INELRYISANSPMLSFYADFSLVNFALPEYPETQTWAEEATLKSLKYRPHSATYRIALYL
		MRQGKVAEAKQWMRATQSYYPYLMPRYADEIRKLPVWAPLLPELLKDCKAFAAAPGHPEA
		KPCK-'

		@expString = 'LIGS+PASTL+LLQSCI'
		@expectedMutTemplate = "atgcccgctgaaacgaccgtatccggcgcgcacccc
		gccgccaaactgccgatttacatcctgccctgcttcctttggataggcatcgtccccttt
		accttcgcgctcaaactgaaaccgtcgcccgacttttaccacgatgccgccgccgcagcc
		ggcctgattgtcctgttgttcctcacggcaggaaaaaaactgtttgatgtcaaaatcccc
		gccatcagcttccttctgtttgcaatggcggcgttttggtatcttcaggcacgcctgatg
		aacctgatttaccccggtatgaacgacatcgtctcttggattttcatcttgctcgccgtc
		agcgcgtgggcctgccggagcttggtcgcacacttcggacaagaacgcatcgtgaccctg
		tttgcctggtcgctgcttatcggctccccggcgagcaccctgctgcttcaatcctgcatc
		gtcgtcatccagtttgccggctgggaagacacccctctgtttcaaaacatcatcgtttac
		agcgggcaaggcgtaatcggacacatcgggcagcgcaacaacctcggacactacctcatg
		tggggcatactcgccgccgcctacctcaacggacaacgaaaaatccccgccgccctcggc
		gtaatctgcctgattatgcagaccgccgttttaggtttggtcaactcgcgcaccatcttg
		acctacatagccgccatcgccctcatccttcccttctggtatttccgttcggacaaatcc
		aacaggcggacgatgctcggcatagccgcagccgtattccttaccgcgctgttccaattt
		tccatgaacaccattctggaaacctttactggcatccgctacgaaactgccgtcgaacgc
		gtcgccaacggcggtttcacagacttgccgcgccaaatcgaatggaataaagcccttgcc
		gccttccagtccgccccgatattcgggcacggctggaacagttttgcccaacaaaccttc
		ctcatcaatgccgaacagcacaacatatacgacaacctcctcagcaacttgttcacccat
		tcccacaacatcgtcctccaactccttgcagagatgggaatcagcggcacgcttctggtt
		gccgcaaccctgctgacgggcattgccgggctgcttaaacgccccctgacccccgcatcg
		cttttcctaatctgcacgcttgccgtcagtatgtgccacagtatgctcgaatatcctttg
		tggtatgtctatttcctcatccctttcggactgatgctcttcctgtcccccgcagaggct
		tcagacggcatcgccttcaaaaaagccgccaatctcggcatactgaccgcctccgccgcc
		atattcgcaggattgctgcacttggactggacatacacccggctggttaacgccttttcc
		cccgccactgacgacagtgccaaaaccctcaaccggaaaatcaacgagttgcgctatatt
		tccgcaaacagtccgatgctgtccttttatgccgacttctccctcgtaaacttcgccctg
		ccggaataccccgaaacccagacttgggcggaagaagcaaccctcaaatcactaaaatac
		cgcccccactccgccacctaccgcatcgccctctacctgatgcggcaaggcaaagttgca
		gaagcaaaacaatggatgcgggcgacacagtcctattacccctacctgatgccccgatac
		gccgacgaaatccgcaaactgcccgtatgggcgccgctgctacccgaactgctcaaagac
		tgcaaagccttcgccgccgcgcccggtcatccggaagcaaaaccctgcaaatga".strip.delete("\n\t")
		@ExperimentObject = InsertionExperiment.new(@expString, @template)
	end

	def test_mutatedTemplate
		assert_not_equal(@template, @ExperimentObject.mutatedTemplate)
	end

	def test_primers
		assert(@ExperimentObject.mutatedTemplate.include?(@ExperimentObject.forwardPrimer.to_s))
		revPrimer = Bio::Sequence::NA.new(@ExperimentObject.reversePrimer)
		assert(@ExperimentObject.mutatedTemplate.include?(revPrimer.reverse_complement))
	end

	def test_insertions
		assert(! @ExperimentObject.template.to_s.include?(@ExperimentObject.insertionSeq))
		assert(@ExperimentObject.mutatedTemplate.include?(@ExperimentObject.insertionSeq))
	end
end





