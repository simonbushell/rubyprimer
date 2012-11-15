require "./Snippet"
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
end
