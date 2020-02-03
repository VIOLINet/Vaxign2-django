function openHTab(evt, tabName) {
	var i, tabcontent, tablinks;
	if ( document.getElementById(tabName).style.display == "block" ) {
		closeHTab(evt, tabName);
	} else {
		tabcontent = document.getElementsByClassName("htabcontent");
		for (i = 0; i < tabcontent.length; i++) {
		    tabcontent[i].style.display = "none";
		}
		tablinks = document.getElementsByClassName("htablinks");
		for (i = 0; i < tablinks.length; i++) {
		    tablinks[i].className = tablinks[i].className.replace(" active", "");
		}
		document.getElementById(tabName).style.display = "block";
		evt.currentTarget.className += " active";
	}
}

function closeHTab(evt, tabName) {
	document.getElementById(tabName).style.display='none';
	tablinks = document.getElementsByClassName("htablinks");
	for (i = 0; i < tablinks.length; i++) {
	    tablinks[i].className = tablinks[i].className.replace(" active", "");
	}
}

function openVTab(evt, tabName) {
	var i, tabcontent, tablinks;
	if ( document.getElementById(tabName).style.display == "block" ) {
		closeVTab(evt, tabName);
	} else {
		tabcontent = document.getElementsByClassName("vtabcontent");
		for (i = 0; i < tabcontent.length; i++) {
		    tabcontent[i].style.display = "none";
		}
		tablinks = document.getElementsByClassName("vtablinks");
		for (i = 0; i < tablinks.length; i++) {
		    tablinks[i].className = tablinks[i].className.replace(" active", "");
		}
		document.getElementById(tabName).style.display = "block";
		evt.currentTarget.className += " active";
	}
}

function closeVTab(evt, tabName) {
	document.getElementById(tabName).style.display='none';
	tablinks = document.getElementsByClassName("vtablinks");
	for (i = 0; i < tablinks.length; i++) {
	    tablinks[i].className = tablinks[i].className.replace(" active", "");
	}
}

