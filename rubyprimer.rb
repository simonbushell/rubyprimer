require './lib/Snippet'
require './lib/Experiments'
require 'sinatra'
require 'sinatra/reloader'
require 'rack'

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
end



