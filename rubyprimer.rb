require './lib/Snippet'
require './lib/Experiments'
require 'sinatra'
#require 'sinatra/reloader'
require 'rack'
require 'bio'
require 'logger'

class RubyPrimerApp < Sinatra::Base

	@@AAcodes = 'acdefghiklmnpqrstvwy'
	
	set :environment, :production

	configure :production, :development do
    	enable :logging, :sessions
    	$logger = Logger.new(STDOUT)
  	end
	
	# configure :development do
 #    	register Sinatra::Reloader
 #  	end

	get '/' do 
		erb :index
	end

	get '/about' do
		erb :about
	end

	post '/submit' do
		experiments = params[:experimentStrings].gsub(/\s+/, "").split(',')
		params[:DNAinput] = params[:DNAinput].gsub(/[^actg]/i, '').downcase
		@results = []
		# $logger.info("input received: #{params[:DNAinput]}\n#{experiments}")
		experiments.each do |e|
			begin
				if e =~ /\A[#{@@AAcodes}]+-[#{@@AAcodes}]+\z/i
					@results << DeletionExperiment.new(e, params[:DNAinput])
				elsif e =~ /\A[#{@@AAcodes}]+\+[#{@@AAcodes}]+\+[#{@@AAcodes}]+\z/i
					@results << InsertionExperiment.new(e, params[:DNAinput])
				elsif e =~ /\A[#{@@AAcodes}]+\*[#{@@AAcodes}]+\*[#{@@AAcodes}]+\z/i
					@results << SubstitutionExperiment.new(e, params[:DNAinput])
				else
					@results << ErrorExperiment.new(e, params[:DNAinput], errorString="Unknown Experiment String")
					$logger.error("#{@results[-1].log}")

				end
			rescue SnippetError, ExperimentError, DNAIndexError
				$logger.error("#{$!}. ExperimentString: #{e}. Template: #{params[:DNAinput]}")
				@results << ErrorExperiment.new(e, params[:DNAinput], errorString=$!)
			rescue RubyPrimerError
				$logger.fatal("RubyPrimerError!! ExperimentString: #{e}. Template: #{params[:DNAinput]}")
				@results << ErrorExperiment.new(e, params[:DNAinput], errorString=$!)
			end
			if @results[-1].class != ErrorExperiment
				$logger.info("Data Submitted: #{@results[-1].log}")
			end
		end
		erb :output 
	end

	post '/ajaxtranslate' do
		dnaobj = Bio::Sequence::NA.new(params[:DNAinput].gsub(/[^actg]/i, '').downcase)
		return dnaobj.translate
	end
end