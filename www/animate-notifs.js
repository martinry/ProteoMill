
Shiny.addCustomMessageHandler('set-icon', function(icon) {
    
  var class_name = "";
  
  if (icon == "success") {
    class_name = "fa fa-check-circle text-success";
  } else if (icon == "info") {
    class_name = "fa fa-info-circle text-info";
  } else if (icon == "danger") {
    class_name = "fa fa-warning text-danger";
  } else {return false}
  
  $('#notif-icon').attr('class',class_name);


  anime.timeline({loop: false})
    .add({
      targets: '#notif-icon',
      opacity: 1,
      easing: "easeOutExpo",
      duration: 1500
    }).add({
      targets: '#notif-icon',
      opacity: 0,
      easing: "easeOutExpo",
      duration: 1500,
      delay: 1400,
    });
  
});


Shiny.addCustomMessageHandler('background-color', function(message) {
  document.querySelector('.letters2').innerText = message;
  document.querySelector('.ml11').style = 'opacity: 1;';

  // Wrap every letter in a span
  $('.ml11 .letters2').each(function(){
    $(this).html($(this).text().replace(/([^\x00-\x80]|\w)/g, "<span class='letter2'>$&</span>"));
  });
  
  var basicTimeline = anime.timeline({
  loop: false
  });
  
//  anime.timeline({loop: false})

basicTimeline
    .add({
      targets: '.ml11 .line',
      scaleY: [0,1],
      opacity: [0.5,1],
      easing: "easeOutExpo",
      duration: 700
    })
    .add({
      targets: '.ml11 .line',
      translateX: [0,$(".ml11 .letters2").width()],
      easing: "easeOutExpo",
      duration: 700,
      delay: 100
    }).add({
      targets: '.ml11 .letter2',
      opacity: [0,1],
      easing: "easeOutExpo",
      duration: 600,
      offset: '-=775',
      delay: function(el, i) {
        return 34 * (i+1);
      }
    }).add({
      targets: '.ml11',
      opacity: 0,
      duration: 900,
      easing: "easeOutExpo",
      delay: 900
    });
    
    
    
//.add({
//      targets: '.letters2',
//      delay: 1000,
//      update: function(){
//        $('.letter2').innerText = '';
//      }
//      
//    });

//  $('.letter2').innerText = "";
//  $('.letter2').empty();
        
      
    
});