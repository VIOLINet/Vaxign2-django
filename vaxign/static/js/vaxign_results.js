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

function loadPopCov(queryID, seqID) {
    if (document.getElementById('population_coverage').innerHTML=='') {
        document.getElementById('population_coverage').innerHTML = "<strong>Opening this page for the first time may take a while.</strong><div class='loading' style='position:relative;'></div>";
        dojo.xhrGet({
            url: "/vaxign2/query/"+queryID+"/protein/"+seqID+"/population_coverage/combine",
            load: function(data){
                $('#population_coverage').html(data);
            },
            error: function(data) {
                alert("An error occurred: " + data);
            },
        });
    };
}

function loadOrthoMCLPhylogeny(queryID, seqID) {
    if (document.getElementById('orthomcl_phylogeny').innerHTML=='') {
        document.getElementById('orthomcl_phylogeny').innerHTML = "<strong>Opening this page for the first time may take a while.</strong><div class='loading' style='position:relative;'></div>";
        dojo.xhrGet({
            url: "/vaxign2/query/"+queryID+"/protein/"+seqID+"/orthomcl/phylogeny",
            load: function(data){
                $('#orthomcl_phylogeny').html(data);
            },
            error: function(data) {
                alert("An error occurred: " + data);
            },
        });
    };
}