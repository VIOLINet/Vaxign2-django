var mhc_i_supertype_alleles = ['HLA-A*01:01', 'HLA-A*30:01', 'HLA-A*02:01', 'HLA-A*03:01', 'HLA-A*24:02', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*44:02', 'HLA-B*58:01', 'HLA-B*15:01'];

var mhc_ii_supertype_alleles = ['HLA-DRB1*01:01','HLA-DPA1*01:03/DPB1*02:01', 'HLA-DQA1*05:01/DQB1*03:01', 'HLA-DRB1*04:01', 'HLA-DRB1*03:01', 'HLA-DPA1*02:01/DPB1*01:01', 'HLA-DQA1*05:01/DQB1*02:01', 'HLA-DRB1*07:01'];

var mhc_i_reference_alleles = ['HLA-A*01:01', 'HLA-A*02:01', 'HLA-A*02:03', 'HLA-A*02:06', 'HLA-A*03:01', 'HLA-A*11:01', 'HLA-A*23:01', 'HLA-A*24:02', 'HLA-A*26:01', 'HLA-A*30:01', 'HLA-A*30:02', 'HLA-A*31:01', 'HLA-A*32:01', 'HLA-A*33:01', 'HLA-A*68:01', 'HLA-A*68:02', 'HLA-B*07:02', 'HLA-B*08:01', 'HLA-B*15:01', 'HLA-B*35:01', 'HLA-B*40:01', 'HLA-B*44:02', 'HLA-B*44:03', 'HLA-B*51:01', 'HLA-B*53:01', 'HLA-B*57:01', 'HLA-B*58:01'];

var mhc_ii_reference_alleles = ['HLA-DPA1*01:03/DPB1*02:01', 'HLA-DPA1*01:03/DPB1*04:01', 'HLA-DPA1*02:01/DPB1*01:01', 'HLA-DPA1*02:01/DPB1*05:01', 'HLA-DPA1*03:01/DPB1*04:02', 'HLA-DQA1*01:01/DQB1*05:01', 'HLA-DQA1*01:02/DQB1*06:02', 'HLA-DQA1*05:01/DQB1*02:01', 'HLA-DQA1*05:01/DQB1*03:01', 'HLA-DRB1*01:01', 'HLA-DRB1*03:01', 'HLA-DRB1*04:01', 'HLA-DRB1*04:05', 'HLA-DRB1*07:01', 'HLA-DRB1*08:02', 'HLA-DRB1*09:01', 'HLA-DRB1*11:01', 'HLA-DRB1*12:01', 'HLA-DRB1*13:02', 'HLA-DRB1*15:01', 'HLA-DRB3*01:01', 'HLA-DRB3*02:02', 'HLA-DRB4*01:01', 'HLA-DRB5*01:01'];

function change_option_host() {
    var mhc_species = document.getElementById('mhc_species');
    var mhc_class = document.getElementById('mhc_class');
    var mhc_allele = document.getElementById('mhc_allele');
    var epitope_length = document.getElementById('epitope_length');
    
    species = mhc_species.options[mhc_species.selectedIndex].value;
    mhc = mhc_class.options[mhc_class.selectedIndex].value;
    
    if (species=='rhesus macaque') {
    	mhc_class.length=1;
    } else {
    	mhc_class.length=2;
    	mhc_class.options[1].value = 'II';
    	mhc_class.options[1].text = 'MHC-II'
    }
    
    mhc_class.disabled = true;
    mhc_allele.disabled = true;
    epitope_length.disabled = true;
    
    dojo.xhrGet({
        url: "/vaxign2/api/t_vaxign_allele_group/"+species+'/'+mhc,
        handleAs: "json",
        load: function(data) {
            var alleleSet = new Set();
            for(i=0; i<data.length; i++) {
                alleleSet.add(data[i].mhc_allele);
            };
            alleleList = Array.from(alleleSet);
            mhc_allele.length = alleleList.length+1;
            mhc_allele.options[0].value = 'any';
            mhc_allele.options[0].text = 'Any Allele';
            for(i=1; i<alleleList.length+1; i++) {
            	mhc_allele.options[i].value = alleleList[i-1].replace(/[^A-Za-z0-9]/g, "");
            	mhc_allele.options[i].text = alleleList[i-1];
            }
        },
        error: function(data) {
            alert("An error occurred: " + data);
        },
        timeout:5000,
    });
    epitope_length.length = 1;
    epitope_length.options[0].value = 'any';
    epitope_length.options[0].text = 'Any Length';
    
    mhc_allele.value = 'any';
    epitope_length.value = 'any';
    
    mhc_class.disabled = false;
    mhc_allele.disabled = false;
    epitope_length.disabled = false;
}

function change_option_class() {
    var mhc_species = document.getElementById('mhc_species');
    var mhc_class = document.getElementById('mhc_class');
    var mhc_allele = document.getElementById('mhc_allele');
    var epitope_length = document.getElementById('epitope_length');
    
    species = mhc_species.options[mhc_species.selectedIndex].value;
    mhc = mhc_class.options[mhc_class.selectedIndex].value;
    
    mhc_species.disabled = true;
    mhc_allele.disabled = true;
    epitope_length.disabled = true;
    
    dojo.xhrGet({
        url: "/vaxign2/api/t_vaxign_allele_group/"+species+'/'+mhc,
        handleAs: "json",
        load: function(data) {
            var alleleSet = new Set();
            for(i=0; i<data.length; i++) {
                alleleSet.add(data[i].mhc_allele);
            };
            alleleList = Array.from(alleleSet);
            mhc_allele.length = alleleList.length+1;
            mhc_allele.options[0].value = 'any';
            mhc_allele.options[0].text = 'Any Allele';
            for(i=1; i<alleleList.length+1; i++) {
            	mhc_allele.options[i].value = alleleList[i-1].replace(/[^A-Za-z0-9]/g, "");
            	mhc_allele.options[i].text = alleleList[i-1]
            }
        },
        error: function(data) {
            alert("An error occurred: " + data);
        },
        timeout:5000,
    });
    epitope_length.length = 1;
    epitope_length.options[0].value = 'any';
    epitope_length.options[0].text = 'Any Length';
    
    mhc_allele.value = 'any';
    epitope_length.value = 'any';
    
    mhc_species.disabled = false;
    mhc_allele.disabled = false;
    epitope_length.disabled = false;
}

function change_option_allele() {
    var mhc_species = document.getElementById('mhc_species');
    var mhc_class = document.getElementById('mhc_class');
    var mhc_allele = document.getElementById('mhc_allele');
    var epitope_length = document.getElementById('epitope_length');
    
    species = mhc_species.options[mhc_species.selectedIndex].value;
    mhc = mhc_class.options[mhc_class.selectedIndex].value;
    
    mhc_species.disabled = true;
    mhc_class.disabled = true;
    epitope_length.disabled = true;
    
    dojo.xhrGet({
        url: "/vaxign2/api/t_vaxign_allele_group/"+species+'/'+mhc,
        handleAs: "json",
        load: function(data) {
            var lengthSet = new Set();
            for(i=0; i<data.length; i++) {
                if (data[i].mhc_allele==mhc_allele.options[mhc_allele.selectedIndex].text) lengthSet.add(data[i].epitope_length);
            };
            lengthList = Array.from(lengthSet);
            epitope_length.length = lengthList.length+1;
            epitope_length.options[0].value = 'any';
            epitope_length.options[0].text = 'Any Length';
            for(i=1; i<lengthList.length+1; i++) {
            	epitope_length.options[i].value = lengthList[i-1]
            	epitope_length.options[i].text = lengthList[i-1]
            }
        },
        error: function(data) {
            alert("An error occurred: " + data);
        },
        timeout:5000,
    });
    
    epitope_length.value = 'any';
    
    mhc_species.disabled = false;
    mhc_class.disabled = false;
    epitope_length.disabled = false;
}

function add_allele() {
    var mhc_species = document.getElementById('mhc_species');
    var mhc_class = document.getElementById('mhc_class');
    var mhc_allele = document.getElementById('mhc_allele');
    var epitope_length = document.getElementById('epitope_length');
    
    var mhc_species_text = mhc_species.options[mhc_species.selectedIndex].text;
    var mhc_class_text = mhc_class.options[mhc_class.selectedIndex].text;
    var mhc_allele_text = mhc_allele.options[mhc_allele.selectedIndex].text;
    var epitope_length_text = epitope_length.options[epitope_length.selectedIndex].text;
    
    var mhc_species_value = mhc_species.options[mhc_species.selectedIndex].value;
    var mhc_class_value = mhc_class.options[mhc_class.selectedIndex].value;
    var mhc_allele_value = mhc_allele.options[mhc_allele.selectedIndex].value;
    var epitope_length_value = epitope_length.options[epitope_length.selectedIndex].value;
	
	var container = document.getElementById('selection_container');
	
	if (mhc_allele_value == 'any') {
		dojo.xhrGet({
	        url: "/vaxign2/api/t_vaxign_allele_group/"+mhc_species_value+'/'+mhc_class_value,
	        handleAs: "json",
	        load: function(data) {
	            var lengthSet = new Set();
	            for(i=0; i<data.length; i++) {
	            	var allele_format = data[i].mhc_allele;
	            	allele_format = allele_format.replace(/[^A-Za-z0-9]/g, "");
	            	var selection1 = document.getElementById(mhc_species_value+'_'+mhc_class_value+'_'+allele_format+'_'+data[i].epitope_length);
	                if (selection1) 
	                	selection1.remove();
	                var selection2 = document.getElementById(mhc_species_value+'_'+mhc_class_value+'_'+allele_format+'_any');
	                if (selection2)
	                	selection2.remove();
	            };
	        },
	        error: function(data) {
	            alert("An error occurred: " + data);
	        },
	        timeout:5000,
	    });
	}
	if (epitope_length_value == 'any') {
		dojo.xhrGet({
	        url: "/vaxign2/api/t_vaxign_allele_group/"+mhc_species_value+'/'+mhc_class_value,
	        handleAs: "json",
	        load: function(data) {
	            var lengthSet = new Set();
	            for(i=0; i<data.length; i++) {
	            	var selection = document.getElementById(mhc_species_value+'_'+mhc_class_value+'_'+mhc_allele_value+'_'+data[i].epitope_length);
	                if (selection && data[i].mhc_allele==mhc_allele_text) 
	                	selection.remove();
	            };
	        },
	        error: function(data) {
	            alert("An error occurred: " + data);
	        },
	        timeout:5000,
	    });
	}

	epitope_any_container = document.getElementById(mhc_species_value+'_'+mhc_class_value+'_'+mhc_allele_value+'_any');
	allele_any_container = document.getElementById(mhc_species_value+'_'+mhc_class_value+'_any_any');
	if (allele_any_container) {
		var message = allele_any_container.textContent;
		alert(message.substring(0,message.length-1)+' is currently selected.');
		return true;
	}
	if (epitope_any_container) { 
		var message = epitope_any_container.textContent;
		alert(message.substring(0,message.length-1)+' is currently selected.');
		return true;
	}
	container_id = mhc_species_value+'_'+mhc_class_value+'_'+mhc_allele_value+'_'+epitope_length_value;
	container_text = mhc_species_text+' | '+mhc_class_text+' | '+mhc_allele_text+' | '+epitope_length_text;
	if (container.innerHTML == '') container.innerHTML = "<span onclick=\"document.getElementById('selection_container').innerHTML='';\" class='topright'>&times</span>";
	if (document.getElementById(container_id) == null)
		container.innerHTML += '<span id="'+container_id+'" class="selected_allele">'+container_text+' <button class="btn">&times</button></span>';
	
	return true;
}

function refresh_table(query_id, sequence_id) {
	if (sequence_id == '') {
		var selection_container = document.getElementById('selection_container');
		
		var groupSet = new Set();
		selection_container.childNodes.forEach( function(selection) {
			if (typeof(selection.id) != "undefined")
				groupSet.add(selection.id);
		});
		
		var groups = '';
		groupSet.forEach( function(group) {
			groups += group + ','
		});
		
		if (groupSet.size > 10 || groupSet.has('human_I_any_any') || groupSet.has('human_II_any_any') )
			if (!confirm("Currently selected over 10(or 'any') alleles . This request may take a while. Do you want to proceed?"))
				return;

		$("#loading").show();
		
		dojo.xhrGet({
	        url: "/vaxign2/api/t_vaxign_mast_results/"+query_id,
	        handleAs: "json",
	        content: {groups: groups},
	        load: function(data) {
	        	var table = $("#vaxitop_results").DataTable()
	        	table.clear();
	        	for(i=0; i<data.length; i++) {
	        		table.row.add([
		               i+1,
		               data[i].protein,
		               data[i].epitope,
		               data[i].mhc_allele,
		               data[i].c_hit_p_value,
		               data[i].c_hit_end-data[i].c_hit_start,
		               data[i].c_hit_start,
		               data[i].c_hit_end,
	        		]).draw();
	        	};
	        	$('td.dataTables_empty').html('No epitopes found.');
	        	$("#loading").hide();
	        },
	        error: function(data) {
	            alert("An error occurred: " + data);
	        },
	    });
	} else {
		var selection_container = document.getElementById('selection_container');
		
		var groupSet = new Set();
		selection_container.childNodes.forEach( function(selection) {
			if (typeof(selection.id) != "undefined")
				groupSet.add(selection.id);
		});
		
		var groups = '';
		groupSet.forEach( function(group) {
			groups += group + ','
		});
		
		if (groupSet.size > 20 || groupSet.has('human_I_any_any') || groupSet.has('human_II_any_any') )
			if (!confirm("Currently selected over 20(or 'any') alleles. This request may take a while. Do you want to proceed?"))
				return;
		
		$("#loading").show();
		
		dojo.xhrGet({
	        url: "/vaxign2/api/t_vaxign_mast_results/"+query_id+'/'+sequence_id,
	        handleAs: "json",
	        content: {groups: groups},
	        load: function(data) {
	        	
	        	
	        	$("#vaxitop_results").DataTable().destroy();
	        	
	        	$.fn.dataTable.ext.search.push(
        		    function( settings, data, dataIndex ) {
        		    	if ( settings.nTable.id !== 'vaxitop_results' ) {
        		            return true;
        		        }
        		        var max = parseFloat( $('#p_val_cut').val() );
        		        var pval = parseFloat( data[3] ) || 0;
        		 
        		        if ( isNaN( max ) || pval <= max )
        		        {
        		            return true;
        		        }
        		        return false;
        		    }
        		);
	        	
	            var table = $('#vaxitop_results').DataTable( {
	                "lengthMenu":[25, 50, 100, 250, 500, 750, 1000],
	                buttons: [
	        	        'copy',
	        	        'csv',
	        	        'excel',
	        	        {extend : 'pdfHtml5',
	        	            title : function() {
	        	                return "Vaxign Vaccine Design\nVaxitop MHC-I and MHC-II Epitope Prediction";
	        	            },
	        	            pageSize : 'LEGAL',
	        	            text : 'PDF',
	        	            titleAttr : 'PDF'
	                } ],
	                "columns": [
	                    {"data": "index"},
	                    {"data": "epitope"},
	                    {"data": "allele"},
	                    {"data": "p_value"},
	                    {"data": "length"},
	                    {"data": "start"},
	                    {"data": "end"},
	                    {
	                         "class":          "details-control",
	                         "orderable":      false,
	                         "data":           null,
	                         "defaultContent": ""
	                     },
	                ],
	                select: {
	                    style: 'multi',
	                    selector: 'td:not(:last-child)',
	                },
	                fixedHeader: true,
	                "dom": 'f<"export vaxitop_export"B><"clear">tli<"vaxitop_page"p>',
	                "language": {
	                    "emptyTable": "No epitopes found."
	                },
	                "fnPreDrawCallback":function(){
	                    $("#loading").show();
	                },
	                "fnDrawCallback":function(){
	                    $("#loading").hide();
	                },
	            } );
			    
	        	for(i=0; i<data.length; i++) {
	        		table.row.add({
	 	               "index": i+1,
	 	               "epitope": data[i].epitope,
	 	               "allele": data[i].mhc_allele,
	 	               "p_value": data[i].c_hit_p_value,
	 	               "length": data[i].c_hit_end-data[i].c_hit_start,
	 	               "start": data[i].c_hit_start,
	 	               "end": data[i].c_hit_end,
	         		}).draw();
	        	};
	        	
	            var detailRows = [];
	            
	            $('#vaxitop_results tbody').on( 'click', 'tr td.details-control', function () {
	                var tr = $(this).closest('tr');
	                var row = table.row( tr );
	                var idx = $.inArray( tr.attr('id'), detailRows );
	         
	                if ( row.child.isShown() ) {
	                    tr.removeClass( 'details' );
	                    row.child.hide();
	         
	                    // Remove from the 'open' array
	                    detailRows.splice( idx, 1 );
	                }
	                else {
	                    tr.addClass( 'details' );
	                    row.child( format( row.data() ) ).show();
	         
	                    // Add to the 'open' array
	                    if ( idx === -1 ) {
	                        detailRows.push( tr.attr('id') );
	                    }
	                }
	            } );
	            
	            table.on( 'draw', function () {
	                $.each( detailRows, function ( i, id ) {
	                    $('#'+id+' td.details-control').trigger( 'click' );
	                } );
	            } );
	            
	            $('#p_val_cut').keyup( function() {
	                table.draw();
	            } );
	            
	        	$("#loading").hide();
	        },
	        error: function(data) {
	            alert("An error occurred: " + data);
	        },
	    });
	};
}

function add_option(option) {
	var container = document.getElementById('selection_container');
	
	if (option=='Human MHC-I Supertype Alleles') {
		container.innerHTML = "<span onclick=\"document.getElementById('selection_container').innerHTML='';\" class='topright'>&times</span>";;
		mhc_i_supertype_alleles.forEach( function(allele) {
			container_id = 'human_I_'+allele.replace(/[^A-Za-z0-9]/g, "")+'_any';
			container_text = 'Human | MHC-I | '+allele+' | Any';
			container.innerHTML += '<span id="'+container_id+'" class="selected_allele">'+container_text+' <button class="btn">&times</button></span>';
		} );
	} else if (option=='Human MHC-II Supertype Alleles') {
		container.innerHTML = "<span onclick=\"document.getElementById('selection_container').innerHTML='';\" class='topright'>&times</span>";;
		mhc_ii_supertype_alleles.forEach( function(allele) {
			container_id = 'human_II_'+allele.replace(/[^A-Za-z0-9]/g, "")+'_any';
			container_text = 'Human | MHC-II | '+allele+' | Any';
			container.innerHTML += '<span id="'+container_id+'" class="selected_allele">'+container_text+' <button class="btn">&times</button></span>';
		} );
	} else if (option=='Human MHC-I & MHC-II Supertype Alleles') {
		container.innerHTML = "<span onclick=\"document.getElementById('selection_container').innerHTML='';\" class='topright'>&times</span>";;
		mhc_i_supertype_alleles.forEach( function(allele) {
			container_id = 'human_I_'+allele.replace(/[^A-Za-z0-9]/g, "")+'_any';
			container_text = 'Human | MHC-I | '+allele+' | Any';
			container.innerHTML += '<span id="'+container_id+'" class="selected_allele">'+container_text+' <button class="btn">&times</button></span>';
		} );
		mhc_ii_supertype_alleles.forEach( function(allele) {
			container_id = 'human_II_'+allele.replace(/[^A-Za-z0-9]/g, "")+'_any';
			container_text = 'Human | MHC-II | '+allele+' | Any';
			container.innerHTML += '<span id="'+container_id+'" class="selected_allele">'+container_text+' <button class="btn">&times</button></span>';
		} );
	
	} else if (option=='Human MHC-I Reference Alleles') {
		container.innerHTML = "<span onclick=\"document.getElementById('selection_container').innerHTML='';\" class='topright'>&times</span>";;
		mhc_i_reference_alleles.forEach( function(allele) {
			container_id = 'human_I_'+allele.replace(/[^A-Za-z0-9]/g, "")+'_any';
			container_text = 'Human | MHC-I | '+allele+' | Any';
			container.innerHTML += '<span id="'+container_id+'" class="selected_allele">'+container_text+' <button class="btn">&times</button></span>';
		} );
	} else if (option=='Human MHC-II Reference Alleles') {
		container.innerHTML = "<span onclick=\"document.getElementById('selection_container').innerHTML='';\" class='topright'>&times</span>";;
		mhc_ii_reference_alleles.forEach( function(allele) {
			container_id = 'human_II_'+allele.replace(/[^A-Za-z0-9]/g, "")+'_any';
			container_text = 'Human | MHC-II | '+allele+' | Any';
			container.innerHTML += '<span id="'+container_id+'" class="selected_allele">'+container_text+' <button class="btn">&times</button></span>';
		} );
	} else if (option=='Human MHC-I & MHC-II Reference Alleles') {
		container.innerHTML = "<span onclick=\"document.getElementById('selection_container').innerHTML='';\" class='topright'>&times</span>";;
		mhc_i_reference_alleles.forEach( function(allele) {
			container_id = 'human_I_'+allele.replace(/[^A-Za-z0-9]/g, "")+'_any';
			container_text = 'Human | MHC-I | '+allele+' | Any';
			container.innerHTML += '<span id="'+container_id+'" class="selected_allele">'+container_text+' <button class="btn">&times</button></span>';
		} );
		mhc_ii_reference_alleles.forEach( function(allele) {
			container_id = 'human_II_'+allele.replace(/[^A-Za-z0-9]/g, "")+'_any';
			container_text = 'Human | MHC-II | '+allele+' | Any';
			container.innerHTML += '<span id="'+container_id+'" class="selected_allele">'+container_text+' <button class="btn">&times</button></span>';
		} );
	}
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