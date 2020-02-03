function chk_opt(elementID1, elementID2) {
    var elementSel1 = document.getElementById(elementID1);
    var elementSel2 = document.getElementById(elementID2);
    for(j=0; j<elementSel1.options.length; j++) {
        if (elementSel1.options[j].selected) {
            elementSel2.options[j].selected = false;
        }
    }
}

function show_selected_no(elementID1, elementID2) {
    var elementSel1 = document.getElementById(elementID1);
    var elementSel2 = document.getElementById(elementID2);
    var count=0;
    for(j=0; j<elementSel1.options.length; j++) {
        if (elementSel1.options[j].selected) {
            count++;
        }
    }
    elementSel2.innerHTML=count;
}

function loadIEDB(queryID, seqID) {
	if (document.getElementById('iedb_epitope').innerHTML=='') {
		document.getElementById('iedb_epitope').innerHTML = "<div class='loading' style='position:relative;'></div>";
		dojo.xhrGet({
			url: "/vaxign2/query/"+queryID+"/protein/"+seqID+"/iedb/search",
			load: function(data){
				$('#iedb_epitope').html(data);
			},
			error: function(data) {
				alert("An error occurred: " + data);
			},
		});
	};
}

function loadVaxitop(queryID, seqID) {
	if (document.getElementById('vaxitop').innerHTML=='') {
		document.getElementById('vaxitop').innerHTML = "<div class='loading' style='position:relative;'></div>";
		dojo.xhrGet({
			url: "/vaxign2/query/"+queryID+"/protein/"+seqID+"/vaxitop",
			load: function(data){
				$('#vaxitop').html(data);
			},
			error: function(data) {
				alert("An error occurred: " + data);
			},
		});
	};
}