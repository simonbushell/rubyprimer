##About
RubyPrimer is a webtool for creating mutagenic oligonucleotide primers for use in the Quikchange II Site-Directed Mutagenesis kit. Our group recently published a modified protocol that allows the kit to be used for large-scale insertions and deletions in protein constructs. The protocol also allows for substitutions to be made at a greater efficiency than via the manufacturer's instructions.

Central to the success of modified protocol is a new primer design strategy. Read [Liu et al (2008)](http://www.biomedcentral.com/1472-6750/8/91) for further details. While the strategy is straightforward, it can be time consuming to implement: especially if one requires many constructs to be designed for downstream high-throughput applications. RubyPrimer makes this task easier.
###How to Use
The only input required is the coding DNA for your protein (pasted in the correct frame - click the "Translate" button to check) and an experiment string to represent the mutation, insertion or deletion you are designing. Multiple experiment strings can be submitted: separate each string with a comma (,). All inputs are case-insensitive. If your modification is close to the N or C terminals of your protein, you might want to include some flanking DNA from your backbone plasmid to prevent an error.

###Experiment Strings
####General Considerations
When designing your strings, the string must be unambigious. For instance in the deletion string **FTGY-DERT** the sequences **FTGY** and **DERT** must only appear once in your template. 4 or 5 residues will usually suffice, but pay particular attention if your modofications occur close to poly-histidine tags. To prevent an error, always include the entire poly-His tag in your experiment string. 
####Deletions
Deleted regions of proteins are represented by a hyphen (-). The experiment string **VNAFS-TGYA** will delete all residues between the S and T. 
####Insertions
Insertions strings are represented by a plus sign (+). Thus the experiment string **VNAFS+AAAAAA+HYTF** will insert 6 alanines in between the serine and the histidine. RubyPrimer will automatically generate a codon-optimised DNA sequence to cater to your insertion. 
####Substitutions
Substitution strings are represented by asterisks (\*). The experiment string **VNAFS\*C\*GTHY** will *substitute* whatever residue follows the serine with a cysteine. RubyPrimer will automatically select the best (ie lowest mismatched) codon to perform the mutation. Note that RubyPrimer can only currently substitute *one* amino acid per experiment. This will likely change in the future. 
###Technical
RubyPrimer is written in the Ruby programming language and deployed on the web via [Sinatra](http://www.sinatrarb.com). RubyPrimer is open-source, free software (via an MIT licence). The source code can be freely viewed at [GitHub](https://github.com/simonbushell/rubyprimer). Contrubutions are welcome. 