


// Wrap every letter in a span
$('.ml9 .letters').each(function(){
  $(this).html($(this).text().replace(/([^\x00-\x80]|\w)/g, "<span class='letter'>$&</span>"));
});

var targetElm = document.querySelector('.ml9');

$( window ).ready(function() {
  setTimeout(function() {
    anime.timeline({loop: false})
      .add({
        targets: '.ml9 .letter',
        scale: [0.25, 1],
        easing: "easeInOutQuart",
        duration: 1400,
        elasticity: 500,
        delay: function(el, i) {
          return 55 * (i+1);
        }
      });
  }, 25);
});

$('.begindiv').delay(1650).queue(function (next) {
    $(this).fadeTo( 1500 , 1);
    next();
});

$('.videoWrapper').delay(800).queue(function (next) {
    $(this).fadeTo( 1500 , 0.85);
    next();
});


$('.begindiv').click(function(){
    $('.overlay').fadeTo(400, 0);
    $('.overlay').queue(function() {
      $( this ).remove();
    });
});



$(document).ready(function() {
  var colors = ["#28313b", "#485461"];
  bg1 = "linear-gradient(180deg," + colors[1] + "," + colors[1] + ")";
  $(".gradient").css("background-image", bg1);

  $(".gradient").mousemove(function(event) {
    var w = $(this).width(),
      pct = 360 * (+event.pageX) / w,
      bg2 = "linear-gradient(" + pct + "deg," + colors[0] + "," + colors[1] + ")";
    $(".gradient").css("background-image", bg2);
  });
  
  
  $(".gradient").mouseleave(function(event) {
    $(".gradient").css("background-image", bg1);
  });
  
  
});


//current = document.querySelector('#notifMenu > a > span').innerHTML;

//if (current !== document.querySelector('#notifMenu > a > span').innerHTML) {
//  alert("hello");
//}

//document.querySelector('#notifMenu > a > span').on('change', function(){
//  document.querySelector('#notifMenu > ul > li:nth-child(2) > ul > li:last-child > a > i');
//});

// document.querySelector('#notifMenu > ul > li:nth-child(2) > ul > li:last-child > a > i');

