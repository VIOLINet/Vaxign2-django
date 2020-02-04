function loadEggNOGFunction(queryID, seqID) {
	if (document.getElementById('eggnog_function').innerHTML=='') {
		document.getElementById('eggnog_function').innerHTML = "<div class='loading' style='position:relative;'></div>";
		dojo.xhrGet({
			url: "/vaxign2/query/"+queryID+"/protein/"+seqID+"/eggnog/function",
			load: function(data){
				$('#eggnog_function').html(data);
			},
			error: function(data) {
				alert("An error occurred: " + data);
			},
		});
	};
}