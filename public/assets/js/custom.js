$(document).ready(function(){

  $(".mutatedTemplate").hide();

  //Event handler to populate div with translation from Ajax
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

  //Event handler to clear DNA from input screen
  $("#clearDNAButton").click(function(){
    $("#dnainput").val("");
    $("#DNAtranslation").html("");
  });  

  //Event handler to show mutated template in output screem
  $('.showMutatedTemplate').click(function(){
    $(this).parents('.outputContainer').find('.mutatedTemplate').slideToggle();
    });
});