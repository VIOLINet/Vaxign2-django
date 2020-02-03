function reloadGenomes(queryUrl) {
    var kw = {
            url: queryUrl,
            handleAs: "json",
            load: function(data){
                elementSel = document.getElementById('id_query_id');
                elementSel.disabled = false;
                elementSel.length = data.length;
//              elementSel.size=data.length;
                elementSel.options[0].disabled = false;
                for(j=0; j<data.length; j++) {
                    elementSel.options[j].value = data[j].c_query_id;
                    elementSel.options[j].text = data[j].c_species_name;
                }
                //initial_options1();
                reloadOptions('/vaxign2/api/t_vaxign_query/ortholog/exclude/'+data[0].c_query_id, 'id_have_orthologs', 'id_have_no_orthologs');
                //alert('');
            },
            error: function(data){
                alert("An error occurred: " + data);
            },
            timeout: 5000
    };
    dojo.xhrGet(kw);
}

function reloadOptions(queryUrl, elementID1, elementID2) {
    var kw = {
            url: queryUrl,
            handleAs: "json",
            load: function(data){
                elementSel = document.getElementById(elementID1);
                elementSel.disabled = false;
                elementSel.length = data.length;
                if (data.length<10) elementSel.size=data.length;
                else elementSel.size=10;
                for(j=0; j<data.length; j++) {
                	elementSel.options[j].selected = false;
                    elementSel.options[j].value = data[j].c_query_id;
                    elementSel.options[j].text = data[j].c_species_name;
                }
                elementSel = document.getElementById(elementID2);
                elementSel.disabled = false;
                elementSel.length = data.length;
                if (data.length<10) elementSel.size=data.length;
                else elementSel.size=10;
                for(j=0; j<data.length; j++) {
                	elementSel.options[j].selected = false;
                    elementSel.options[j].value = data[j].c_query_id;
                    elementSel.options[j].text = data[j].c_species_name;
                }
            },
            error: function(data){
                alert("An error occurred: " + data);
            },
            timeout: 5000
    };
    dojo.xhrGet(kw);
}

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

var sequence1 = ">gi|62317454|ref|YP_223307.1| SodC, superoxide dismutase, Cu-Zn [Brucella abortus biovar 1 str. 9-941]\nMKSLFIASTMVLMAFPAFAESTTVKMYEALPTGPGKEVGTVVISEAPGGLHFKVNMEKLTPGYHGFHVHE\nNPSCAPGEKDGKIVPALAAGGHYDPGNTHHHLGPEGDGHMGDLPRLSANADGKVSETVVAPHLKKLAEIK\nQRSLMVHVGGDNYSDKPEPLGGGGARFACGVIE";
var sequence2 = ">gi|47566476|ref|YP_016495.2| protective antigen [Bacillus anthracis str. 'Ames Ancestor']\nMKKRKVLIPLMALSTILVSSTGNLEVIQAEVKQENRLLNESESSSQGLLGYYFSDLNFQAPMVVTSSTTG\nDLSIPSSELENIPSENQYFQSAIWSGFIKVKKSDEYTFATSADNHVTMWVDDQEVINKASNSNKIRLEKG\nRLYQIKIQYQRENPTEKGLDFKLYWTDSQNKKEVISSDNLQLPELKQKSSNSRKKRSTSAGPTVPDRDND\nGIPDSLEVEGYTVDVKNKRTFLSPWISNIHEKKGLTKYKSSPEKWSTASDPYSDFEKVTGRIDKNVSPEA\nRHPLVAAYPIVHVDMENIILSKNEDQSTQNTDSQTRTISKNTSTSRTHTSEVHGNAEVHASFFDIGGSVS\nAGFSNSNSSTVAIDHSLSLAGERTWAETMGLNTADTARLNANIRYVNTGTAPIYNVLPTTSLVLGKNQTL\nATIKAKENQLSQILAPNNYYPSKNLAPIALNAQDDFSSTPITMNYNQFLELEKTKQLRLDTDQVYGNIAT\nYNFENGRVRVDTGSNWSEVLPQIQETTARIIFNGKDLNLVERRIAAVNPSDPLETTKPDMTLKEALKIAF\nGFNEPNGNLQYQGKDITEFDFNFDQQTSQNIKNQLAELNATNIYTVLDKIKLNAKMNILIRDKRFHYDRN\nNIAVGADESVVKEAHREVINSSTEGLLLNIDKDIRKILSGYIVEIEDTEGLKEVINDRYDMLNISSLRQD\nGKTFIDFKKYNDKLPLYISNPNYKVNVYAVTKENTIINPSENGDTSTNGIKKILIFSKKGYEIG";

var sequence3 = "62317454";
var sequence4 = "YP_016495.2";

var format1="protein_fasta";
var format2="protein_fasta";
var format3="protein_gi";
var format4="protein_refseq";
    
function setSequence(seqid, format_value) {
	document.getElementById('id_bacteria_strain').style.display = 'inline-block';
	if (seqid == sequence1) {
		document.getElementById('id_organism').value = 'bacteria';
		document.getElementById('id_bacteria_strain').value = '--negative';
	} else {
		document.getElementById('id_organism').value = 'bacteria';
		document.getElementById('id_bacteria_strain').value = '--positive';
	}
    document.getElementById('id_sequence').value=seqid;
    document.getElementById('id_sequence_type').value=format_value;
}

function setSequenceOnly(seqid, format_value) {
    document.getElementById('id_sequence').value=seqid;
    document.getElementById('id_sequence_type').value=format_value;
}

function selectDefaultAnalysis() {
	selectAllAnalysis(false);
	for (var i=0; i<3 ; i++) {
    	document.getElementById( "id_basic_analysis_"+i ).checked=true;
    }
}


function selectAllAnalysis(checked) {
    for (var i=0; i<6 ; i++) {
    	document.getElementById( "id_basic_analysis_"+i ).checked=checked;
    }
}

function showStatsTable(e) {
	var className = $(e).attr('class')+'_detail';
	
	if ($('.' + className).hasClass("hide")) {
		$('.' + className).removeClass("hide");
	} else {
		$('.' + className).addClass("hide");
	}
	
}

