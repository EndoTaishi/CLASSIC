var variables = []
var variants = []
var groups = []
var variantsVisited = false;

function Variable(id, standardName, longName, shortName, units, group, bareGround) {
	this.id = id;
	this.standardName = standardName;
	this.longName = longName;
	this.shortName = shortName;
	this.units = units;
	this.group = group;
	this.bareGround = bareGround;
}

function Variant(variable, frequency, form, nic) {
	this.variable = variable;
	this.frequency = frequency;
	this.form = form;
	this.nic = nic;
}

function main() {
	sceneSetup();
}

//Configure the scene
function sceneSetup() {
	$( document ).tooltip();
	$('button').button();
	$('input#inputXmlFile').button();
	$('button#loadXml').on('click', loadXml);
	$('button#addVariable').on('click', addVariable);
	$('button#removeVariable').on('click', removeVariable);
	$('button#addGroup').on('click', addGroupFromForm);
	$('#accordion').accordion({
		collapsible : true
	});

	$('#tabs').tabs({
		activate : function(event, ui) {
			if (variantsVisited) {
				saveVariants();
				variantsVisited = false;
			}
			switch (event.currentTarget.id) {
			case 'tab1':

				break;
			case 'tab2':
				addGroupsToForm();
				break;
			case 'tab3':
				generateVariantsEditor();
				break;
			case 'tab4':
				generateVariants();
				break;
			}
		}
	});
}

function addGroupFromForm() {
	var group = $('input#addGroup')[0].value.trim();
	if (group != '') {
		addGroup(group);
		addGroupsToForm();
	} else
		alert('The group name can\'t be empty');
}

function addGroupsToForm() {
	$('#group').empty();

	for (var i = 0; i < groups.length; i++)
		$('#group').append($('<option/>').text(groups[i]));

}

function addGroup(group) {
	if (groups.indexOf(group) < 0)
		groups.push(group);
}

//Loads an xml file into the input text area
function loadXml() {
	variantsVisited = false;
	var fileToLoad = document.getElementById("inputXmlFile").files[0];
	var fileReader = new FileReader();
	fileReader.onload = function(fileLoadedEvent) {
		groups = [];
		variants = [];
		variables = [];

		var textFromFileLoaded = fileLoadedEvent.target.result;
		$("textarea#inputXml").text(textFromFileLoaded);
		var temp = $.parseXML(textFromFileLoaded);
		var $inputXML = $(temp);

		$('variable', $inputXML)
		.each(
				function() {
					var id = $(this).attr('id');
					var standardName = $('standardName', this).text().trim();
					var longName = $('longName', this).text().trim();
					var shortName = $('shortName', this).text().trim();
					var units = $('units', this).text().trim();
					var group = $(this.parentElement.attributes['type'])[0].value;
					var bareGround = $(this.attributes['includeBareGround'])[0].value;
					addGroup(group);

					var variable = new Variable(id, standardName, longName, shortName, units, group, bareGround);
					variables.push(variable);
				})

				$('variant', $inputXML).each(function() {
					var freq = $('timeFrequency', this).text().trim();
					switch (freq) {
					case 'annually':
						freq = 0;
						break;
					case 'monthly':
						freq = 1;
						break;
					case 'daily':
						freq = 2;
						break;
					case '6 hourly':
						freq = 3;
						break;
					case 'half hourly':
						freq = 4;
						break;
					default:
						freq = 'PROBLEM';
					}

					var form = $('outputForm', this).text().trim();
					switch (form) {
					case 'grid':
						form = 0;
						break;
					case 'pft':
						form = 1;
						break;
					case 'tile':
						form = 2;
						break;
					case 'layer':
						form = 3;
						break;
					}

					var nic = $('nameInCode', this).text().trim();
					var id = $(this.parentElement.attributes['id'])[0].value;

					var variant = new Variant(id, freq, form, nic);
					variants.push(variant);
				})
	};
	fileReader.readAsText(fileToLoad, "UTF-8");
}

//Save XML file
function saveTextAsFile() {
	var textToSave = document.getElementById("outputXml").value;
	var textToSaveAsBlob = new Blob([ textToSave ], {
		type : "text/plain"
	});
	var textToSaveAsURL = window.URL.createObjectURL(textToSaveAsBlob);
	var fileNameToSaveAs = "output.xml";

	var downloadLink = document.createElement("a");
	downloadLink.download = fileNameToSaveAs;
	downloadLink.innerHTML = "Download File";
	downloadLink.href = textToSaveAsURL;
	downloadLink.onclick = destroyClickedElement;
	downloadLink.style.display = "none";
	document.body.appendChild(downloadLink);

	downloadLink.click();
}

function destroyClickedElement(event) {
	document.body.removeChild(event.target);
}

//Adds a variable manually
function addVariable() {
	var standardName = $('#standardName')[0].value;
	var longName = $('#longName')[0].value;
	var shortName = $('#shortName')[0].value;
	var units = $('#units')[0].value;
	var group = $('#group')[0].value;
	var bareGround = $('#bareGround')[0].value;
	if ($('input#bareGround')[0].checked) bareGround = 'true';
	else bareGround = 'false';

	if (shortName != '') {
		var found = false;
		for (var i = 0; i < variables.length; i++) {
			var temp = variables[i];
			if (temp.shortName == shortName) {
				found = true;
				alert('The variable ' + shortName + ' already exists!');
			}
		}
		if (!found) {
			var variable = new Variable(variables.length, standardName, longName, shortName,
					units, group, bareGround);
			variables.push(variable);
			alert('Successfully added the ' + shortName + ' variable!');
		}
	} else
		alert('Please supply at least the short name for the variable you intend to add!');
}

//Remove variable
function removeVariable() {
	var shortName = $('input#removeVariable')[0].value;
	if (shortName != '') {
		var found = false;
		for (var i = 0; i < variables.length; i++) {
			var temp = variables[i];
			if (temp.shortName == shortName) {
				variables.splice(i, 1);
				found = true;
				alert('Successfully removed the ' + shortName + ' variable!');
			}
		}
		if (!found)
			alert('Couldn\'t find ' + shortName
					+ ' within the loaded variables!');
	} else
		alert('Please supply the short name for the variable you intend to remove!');
}

//Generates a variant editor for a certain variable
function generateVariantsEditor() {
	$('#accordion').empty();
	variantsVisited = true;

	if (variables.length > 0) {
		// Add variables
		for (var i = variables.length - 1; i >= 0 ; i--) {
			var temp = variables[i];
			$('#accordion').append($('<h3/>').text(temp.standardName)).append(
					generateForm(temp.id, temp.shortName, temp.standardName,
							temp.longName, temp.units, temp.group));
		}

		// Refresh variants based on existing values
		for (var i = 0; i < variants.length; i++) {
			var variant = variants[i];
			var id = variant.variable * 1000 + variant.frequency * 10 + variant.form;
			$('input#' + id + '.checkbox').prop('checked', true);
			$('input#' + id + '.textInput').prop('disabled', false).prop(
					'value', variant.nic);
		}

		// Refresh structure
		$('#accordion').accordion('refresh');
	} else
		alert('Please add some variables first!');
}

//Build the variant form for a certain variable
function generateForm(baseId, name, standardName, longName, units, group) {
	var $variable = $("<div/>").addClass('variable').attr('id', baseId).append(
			$("<h4/>").text('Short Name: ' + name)).append(
					$("<p/>").text('Standard Name: ' + standardName)).append(
							$("<p/>").text('Long Name: ' + longName)).append(
									$("<p/>").text('Units: ' + units)).append(
											$("<p/>").text('Group: ' + group));

	baseId *= 1000;

	for (var i = 0; i < 5; i++)
		for (var j = 0; j < 4; j++) {
			var tempName = '';
			var id = baseId + i * 10 + j;

			switch (i) {
			case 0:
				tempName += 'Annually';
				break;
			case 1:
				tempName += 'Monthly';
				break;
			case 2:
				tempName += 'Daily';
				break;
			case 3:
				tempName += '6 hourly';
				break;
			case 4:
				tempName += 'Half-hourly';
				break;
			default:
				tempName += 'PROBLEM';
			}

			tempName += ' - ';

			switch (j) {
			case 0:
				tempName += 'Grid';
				break;
			case 1:
				tempName += 'Per PFT';
				break;
			case 2:
				tempName += 'Per tile';
				break;
			case 3:
				tempName += 'Per layer';
				break;
			default:
				tempName += 'PROBLEM';
			}

			var $variant = $('<div/>').append($('<input/>', {
				type : 'checkbox',
				id : id,
				class : 'checkbox',
				click : function() {
					if (this.checked) {
						enable(this.id)
					} else
						disable(this.id)
				}
			}));

			$variant.append(tempName);
			$variant.append($('<input/>', {
				type : 'text',
				id : id,
				class : 'textInput',
				placeholder : 'Name in code',
				disabled : true
			}));

			$variable.append($variant);
		}

	return $variable;
}

function enable(id) {
	$('input#' + id + '.textInput')[0].disabled = false;
}

function disable(id) {
	$('input#' + id + '.textInput')[0].disabled = true;
	$('input#' + id + '.textInput')[0].value = '';
}

function saveVariants() {
	variants = [];
	for (var v = 0; v < variables.length; v++) {
		for (var i = 0; i < 5; i++)
			for (var j = 0; j < 4; j++) {
				id = v * 1000 + i * 10 + j;
				if ($('input#' + id + '.checkbox')[0].checked) {
					var freq = (Math.floor(id / 10)) % 10;
					var form = id % 10;
					var nic = $('input#' + id + '.textInput')[0].value;

					var variant = new Variant(v, freq, form, nic);
					variants.push(variant);
				}
			}
	}
}

function generateVariants() {
	$('#outputViewer').empty();
	$('#outputViewer')
	.append(
			$('<textarea rows="20" cols="100" id="outputXml" readonly></textarea>'));

	if (variants.length > 0) {
		// Create empty XML document
		var xml = $.parseXML('<?xml version="1.0"?><variableSet/>');

		var d = new Date();
		var strDate = d.getFullYear() + "/" + (d.getMonth()+1) + "/" + d.getDate();
		$('variableSet', xml).attr('type', 'CLASS').attr('version', '1.0').attr('created', strDate);

		// Append groups
		for (var g = 0; g < groups.length; g++) {
			var $group = $(xml.createElement('group')).attr('type', groups[g]);
			$('variableSet', xml).append($group);
		}

		// Append variables
		for (var v = 0; v < variables.length; v++) {
			var $variable = $(xml.createElement('variable')).attr('id', variables[v].id);
			var $standardName = $(xml.createElement('standardName')).text(variables[v].standardName);
			var $longName = $(xml.createElement('longName')).text(variables[v].longName);
			var $shortName = $(xml.createElement('shortName')).text(variables[v].shortName);
			var $units = $(xml.createElement('units')).text(variables[v].units);

			$variable.append($standardName).append($longName)
			.append($shortName).append($units).attr('includeBareGround', variables[v].bareGround);
			$('group[type=' + variables[v].group + ']', xml).append($variable);
		}

		// Append variants
		for (var i = 0; i < variants.length; i++) {
			var variable = variants[i].variable;
			var freq = variants[i].frequency;
			var form = variants[i].form;
			var nic = variants[i].nic;

			var $freq = $(xml.createElement('timeFrequency'));
			var $form = $(xml.createElement('outputForm'));
			var $nic = $(xml.createElement('nameInCode'));

			switch (freq) {
			case 0:
				$freq.html('annually');
				break;
			case 1:
				$freq.html('monthly');
				break;
			case 2:
				$freq.html('daily');
				break;
			case 3:
				$freq.html('6 hourly');
				break;
			case 4:
				$freq.html('half hourly');
				break;
			}

			switch (form) {
			case 0:
				$form.html('grid');
				break;
			case 1:
				$form.html('pft');
				break;
			case 2:
				$form.html('tile');
				break;
			case 3:
				$form.html('layer');
				break;
			}

			$nic.html(nic);

			var $variant = $(xml.createElement('variant'));
			$variant.append($freq).append($form).append($nic);

			$('variable[id=' + variable + ']', xml).append($variant);
		}

		var s = new XMLSerializer();
		var string = s.serializeToString(xml);

		$('#outputXml').text(string).format({
			method : 'xml'
		});

		saveTextAsFile();
	} else
		alert('You must first add some variants!');
}

//Console.log replacement for ease of writing
function print(text) {
	console.log(text);
}

//Upon loading the document
$(document).ready(function() {
	main();
});
