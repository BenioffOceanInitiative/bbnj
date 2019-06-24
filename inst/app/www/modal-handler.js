// This recieves messages of type "modalmessage" from the server.
Shiny.addCustomMessageHandler("modalmessage",
  function(msg) {
    //alert(JSON.stringify(msg));

    $('#modal').find('iframe')
      .prop('src', function(){ return msg.src });

    $('#modal-title').html( msg.title );

    $('#modal').on('show.bs.modal', function () {
      $('.modal-content').css('height',$( window ).height()*0.9);
      $('.modal-body').css('height','calc(100% - 65px - 55.33px)');
    });

    $('#modal').modal();
  }
);
