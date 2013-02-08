require './lib/Snippet'
require './lib/Experiments'
require 'sinatra'
require 'sinatra/reloader'
require 'rack'
require 'Bio'
require 'logger'

class RubyPrimerApp < Sinatra::Base

	@@AAcodes = 'acdefghiklmnpqrstvwy'
	
	configure :production, :development do
    	enable :logging
    	$logger = Logger.new(STDOUT)
  	end
	
	configure :development do
    	register Sinatra::Reloader
  	end

	get '/' do 
		erb :index
	end

	post '/submit' do
		experiments = params[:experimentStrings].gsub(/\s+/, "").split(',')
		@results = []
		$logger.info("input received: #{params[:DNAinput]}\n#{experiments}")
		experiments.each do |e|
			if e =~ /\A[#{@@AAcodes}]+-[#{@@AAcodes}]+\z/i
				@results << DeletionExperiment.new(e, params[:DNAinput])
			elsif e =~ /\A[#{@@AAcodes}]+\+[#{@@AAcodes}]+\+[#{@@AAcodes}]+\z/i
				@results << InsertionExperiment.new(e, params[:DNAinput])
			elsif e =~ /\A[#{@@AAcodes}]+\*[#{@@AAcodes}]+\*[#{@@AAcodes}]+\z/i
				@results << SubstitutionExperiment.new(e, params[:DNAinput])
			else
				$logger.warn("Error Experiment created")
				@results << ErrorExperiment.new(e, params[:DNAinput])
			end			
		end
		erb :output 
	end

	post '/ajaxtranslate' do
		dnaobj = Bio::Sequence::NA.new(params[:DNAinput])
		return dnaobj.translate
	end
end