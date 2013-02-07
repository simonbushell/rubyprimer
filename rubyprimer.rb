require './lib/Snippet'
require './lib/Experiments'
require 'sinatra'
require 'sinatra/reloader'
require 'rack'
require 'Bio'

class RubyPrimerApp < Sinatra::Base

	configure :development do
    	register Sinatra::Reloader
  	end

	get '/' do 
		erb :index
	end

	post '/submit' do
		"#{params}"
	end

	post '/ajaxtranslate' do
		dnaobj = Bio::Sequence::NA.new(params[:DNAinput])
		return dnaobj.translate
	end
end