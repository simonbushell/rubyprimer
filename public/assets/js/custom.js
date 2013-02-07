$(document).ready(function(){
  $("#translateButton").click(function(){
    var text = $("#dnainput").val()
    $.post("/ajaxtranslate",
    {
      DNAinput: text,
    },
    function(data,status){
      $("#DNAtranslation").html("<p>" + data + "</p>");
    });
  });

  $("#clearDNAButton").click(function(){
    $("#dnainput").val("");
    $("#DNAtranslation").html("");
  });  
});