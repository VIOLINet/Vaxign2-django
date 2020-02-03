function changeQueryPublic(queryID) {
	var flag = false;
	var projectID = document.getElementById('project_id').value
	var url;
	if (document.getElementById('query_public').value=='1') {
		if (document.getElementById('project_public').innerText=='Yes') flag = true;
		else if (confirm("The current project is flagged as private.</br>Do you want to make this query public?")) flag = true;
		url = '/vaxign2/project/' + projectID + '/query/' + queryID + '/public';
	} else url = '/vaxign2/project/' + projectID + '/query/' + queryID + '/private';
	dojo.xhrGet({
		url: url,
		error: function(data) {
			alert("An error occurred: " + data);
		},
		timeout:5000,
	});
}

function addCurator(projectID) {
	var user = document.getElementById('curator').querySelector('input').value;
	if (user != '') {
		var url = "/vaxign2/project/" + projectID + "/curator/" + user + "/add";
		dojo.xhrGet({
			url: url,
			error: function(data) {
				alert("An error occurred: " + data);
			},
			timeout:5000,
		});
		location.reload();
	};
}

function removeCurator(projectID, user) {
	var url = "/vaxign2/project/" + projectID + "/curator/" + user + "/remove";
	dojo.xhrGet({
		url: url,
		error: function(data) {
			alert("An error occurred: " + data);
		},
		timeout:5000,
	});
	location.reload();
}

function removeProject(projectID) {
	if (confirm("Are you sure you want to delete this project?")) {
		dojo.xhrGet({
			url: "/vaxign2/project/"+projectID+"/remove",
			error: function(data) {
				alert("An error occurred: " + data);
			},
			timeout:5000,
		});
		window.open("/vaxign2/project","_self");
	}
}

function removeQuery(projectID, queryID) {
	if (confirm("Are you sure you want to delete this query?")) {
		dojo.xhrGet({
			url: "/vaxign2/project/"+projectID+"/query/"+queryID+"/remove",
			error: function(data) {
				alert("An error occurred: " + data);
			},
			timeout:5000,
		});
		location.reload();
	}
}