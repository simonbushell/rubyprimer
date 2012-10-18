require "./rubyprimer"
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
		aaseq = "MAAFWY"
		d = @snippet.backTranslate(aaseq)
		assert_instance_of(Snippet, d)
		assert_equal(d.translate, aaseq)
		assert_equal(d.start, 240)
		assert_equal(d.end, 258)
	end
end
