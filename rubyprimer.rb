require './lib/Snippet'
require './lib/Experiments'
require 'sinatra'
require 'haml'


get '/' do 
	erb :index
end

post '/' do
	m = params[:namea]
	"Hello #{m}"
end

