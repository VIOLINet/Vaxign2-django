function loadEggNOGFunction(queryID, seqID) {
	if (document.getElementById('eggnog_function').innerHTML=='') {
		document.getElementById('eggnog_function').innerHTML = "<strong>Opening this page for the first time may take a while.</strong><div class='loading' style='position:relative;'></div>";
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

function loadEggNOGOrtholog(queryID, seqID) {
	if (document.getElementById('eggnog_ortholog').innerHTML=='') {
		document.getElementById('eggnog_ortholog').innerHTML = "<strong>Opening this page for the first time may take a while.</strong><div class='loading' style='position:relative;'></div>";
		dojo.xhrGet({
			url: "/vaxign2/query/"+queryID+"/protein/"+seqID+"/eggnog/ortholog",
			load: function(data){
				$('#eggnog_ortholog').html(data);
			},
			error: function(data) {
				alert("An error occurred: " + data);
			},
		});
	};
}

function reloadOrtholog(queryID, seqID) {
	var view = document.getElementById('ortholog_view_option');
	var option = view.options[view.selectedIndex].value;
	document.getElementById('eggnog_ortholog').innerHTML = "<strong>Opening this page for the first time may take a while.</strong><div class='loading' style='position:relative;'></div>";
	dojo.xhrGet({
		url: "/vaxign2/query/"+queryID+"/protein/"+seqID+"/eggnog/ortholog/"+option,
		load: function(data){
			$('#eggnog_ortholog').html(data);
		},
		error: function(data) {
			alert("An error occurred: " + data);
		},
	});
}