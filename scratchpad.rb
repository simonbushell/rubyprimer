require "./Snippet"
require "Bio"
require "amatch"

def getBestCodon(codon, aminoacid)
	include Amatch
	codonTable = Bio::CodonTable[1]
	possibleCodons = codonTable.revtrans(aminoacid)
	codonScores = Hash.new
	m = Sellers.new(codon)
	possibleCodons.each {|x| codonScores[x] = m.match(x)}
	return codonScores.sort_by{|codon, score| score}[0][0]
end
