

// Wrap every letter in a span
$('.ml9 .letters').each(function(){
  $(this).html($(this).text().replace(/([^\x00-\x80]|\w)/g, "<span class='letter'>$&</span>"));
});

anime.timeline({loop: false})
  .add({
    targets: '.ml9 .letter',
    scale: [0, 1],
    duration: 1200,
    elasticity: 600,
    delay: function(el, i) {
      return 45 * (i+1)
    }
  })

$('.begindiv').delay(1650).queue(function (next) {
    $(this).fadeTo( 1500 , 1);
    next();
});


$('.begindiv').click(function(){
    $('.overlay').fadeTo(400, 0);
    $('.overlay').queue(function() {
      $( this ).remove();
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

