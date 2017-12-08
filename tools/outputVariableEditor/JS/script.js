var variables = []
var variants = []
var groups = []
var variablesChanged = false;

function Variable(id, standardName, longName, shortName, units, group, bareGround) {
	this.id = id;
	this.standardName = standardName;
	this.longName = longName;
	this.shortName = shortName;
	this.units = units;
	this.group = group;
	this.bareGround = bareGround;
}

function Variant(variable, frequency, form, nic, units) {
	this.variable = variable;
	this.frequency = frequency;
	this.form = form;
	this.nic = nic;
	this.units = units;
}

function main() {
	sceneSetup();
}

//Configure the scene
function sceneSetup() {
	$( document ).tooltip({track: true});
	$('button').button();
	$('input#inputXmlFile').button();
	$('button#loadXml').on('click', loadXml);
	$('button#addVariable').on('click', addVariable);
	$('button#removeVariable').on('click', removeVariable);
	$('button#addGroup').on('click', addGroupFromForm);
	$( "#clearForm" ).on('click', clearForm)
	$('#accordion').accordion({
		collapsible : true
	});

	$('#tabs').tabs({
		activate : function(event, ui) {
			if (variablesChanged) {
				saveVariableChanges();
				variablesChanged = false;
			}
			switch (event.currentTarget.id) {
			case 'tab1':

				break;
			case 'tab2':
				addGroupsToForm();
				break;
			case 'tab3':
				buildVariableConfigForms();
				//generateVariantsEditor();
				break;
			case 'tab4':
				generateVariants();
				break;
			}
		}
	});
}

function clearForm() {
	$('input#shortName')[0].value = '';
	$('input#standardName')[0].value = '';
	$('input#longName')[0].value = '';
	$('input#units')[0].value = '';
	$('input#bareGround')[0].checked = false;
	
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
	variablesChanged = false;
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
					var bareGround = ($(this.attributes['includeBareGround'])[0].value == 'true');
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
					case 'halfhourly':
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
					//var vunits = $('')
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

	var downloadLink = document.getElementById('downloadButton');
	//downloadLink.download = $('input#filename')[0].value;
	//downloadLink.innerHTML = "Download XML";
	downloadLink.href = textToSaveAsURL;
	//downloadLink.onclick = destroyClickedElement;
	//downloadLink.style.display = "none";
	//document.body.appendChild(downloadLink);
	
	//downloadLink.click();
}

//function destroyClickedElement(event) {
//	document.body.removeChild(event.target);
//}

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

// Builds the variable configuration forms
function buildVariableConfigForms() {
	$('#accordion').empty();
	if (variables.length > 0) {
		// Add variables and variants to the form
		for (var v = variables.length - 1; v >= 0 ; v--) {
			var $line;
			var currentVariable = variables[v];
			
			// Add variable
			var $variableForm = $('<div/>').addClass('variable').attr('id', currentVariable.id);
			
			var $variableTable = $('<table class="variableTable"/>');

			$line = $('<tr/>')
			.append($('<td/>').text('Short Name: '))
			.append($('<td/>').append($('<input/>', {
									type: 'text',
									id : currentVariable.id,
									class : 'variableTextInput shortName',
									value : currentVariable.shortName
								})));
			$variableTable.append($line);
			
			$line = $('<tr/>')
			.append($('<td/>').text('Standard Name: '))
			.append($('<td/>').append($('<input/>', {
									type: 'text',
									id : currentVariable.id,
									class : 'variableTextInput standardName',
									value : currentVariable.standardName
								})));
			$variableTable.append($line);
			
			$line = $('<tr/>')
			.append($('<td/>').text('Long Name: '))
			.append($('<td/>').append($('<input/>', {
									type: 'text',
									id : currentVariable.id,
									class : 'variableTextInput longName',
									value : currentVariable.longName
								})));
			$variableTable.append($line);

			$line = $('<tr/>')
			.append($('<td/>').text('Units Name: '))
			.append($('<td/>').append($('<input/>', {
									type: 'text',
									id : currentVariable.id,
									class : 'variableTextInput units',
									value : currentVariable.units
								})));
			$variableTable.append($line);
			
			$line = $('<tr/>')
			.append($('<td/>').text('Include bare ground: '))
			.append($('<td/>').append($('<input/>', {
									type: 'checkbox',
									id : currentVariable.id,
									class : 'variableTextInput bareGround',
									checked : currentVariable.bareGround
								})));
			$variableTable.append($line);
			
			$variableForm.append($('<div/>').append($variableTable));
			
			// Add variants
			var $variantsTable = $('<table class="variantsTable"/>');
			var baseId = currentVariable.id * 1000;

			for (var i = 0; i < 5; i++)
				for (var j = 0; j < 4; j++) {
					var id = baseId + i * 10 + j;
					var temp = '';

					var $variant = $('<tr/>')
						.append($('<td/>')
								.append($('<input/>', {
									type : 'checkbox',
									id : id,
									class : 'checkbox',
									click : function() {
										if (this.checked) enable(this.id)
										else disable(this.id)
									}
								})));
					
					switch (i) {
					case 0: temp = 'Annually'; break;
					case 1: temp = 'Monthly'; break;
					case 2: temp = 'Daily'; break;
					case 3: temp = '6 hourly'; break;
					case 4: temp = 'Half-hourly'; break;
					default: temp = 'PROBLEM';
					}

					$variant.append($('<td/>').text(temp));

					switch (j) {
					case 0: temp = 'Grid'; break;
					case 1: temp = 'Per PFT'; break;
					case 2: temp = 'Per tile'; break;
					case 3: temp = 'Per layer'; break;
					default: temp = 'PROBLEM';
					}

					$variant.append($('<td/>').text(temp));
					
					$variant.append($('<td/>').append($('<input/>', {
						type : 'text',
						id : id,
						class : 'nameInCode',
						placeholder : 'Name in code',
						disabled : true
					})));
					
					$variant.append($('<td/>').append($('<input/>', {
						type : 'text',
						id : id,
						class : 'variantUnits',
						placeholder : currentVariable.units,
						disabled : true
					})));

					$variantsTable.append($variant);
					$variableForm.append($('<div/>').append($variantsTable));
				}

			$('#accordion')
				.append($('<h3/>').text(currentVariable.shortName + " - " + currentVariable.standardName))
				.append($variableForm);
		}

		// Refresh variants data based on existing values
		for (var i = 0; i < variants.length; i++) {
			var variant = variants[i];
			var id = variant.variable * 1000 + variant.frequency * 10 + variant.form;
			$('input#' + id + '.checkbox').prop('checked', true);
			$('input#' + id + '.nameInCode').prop('disabled', false).prop('value', variant.nic);
			$('input#' + id + '.variantUnits').prop('disabled', false).prop('value', variant.units);
		}
		
		// Refresh accordion structure
		$('#accordion').accordion('refresh');
		variablesChanged = true;
	} else
		alert('Please add some variables first!');
}

function enable(id) {
	$('input#' + id + '.nameInCode')[0].disabled = false;
	$('input#' + id + '.variantUnits')[0].disabled = false;
	//$('input#' + id)[0].disabled = false;
}

function disable(id) {
	$('input#' + id + '.nameInCode')[0].disabled = true;
	$('input#' + id + '.nameInCode')[0].value = '';
	$('input#' + id + '.variantUnits')[0].disabled = true;
	$('input#' + id + '.variantUnits')[0].value = '';
}

function saveVariableChanges() {
	variants = [];
	for (var v = 0; v < variables.length; v++) {
		variables[v].shortName = $('input#' + variables[v].id + '.shortName')[0].value;
		variables[v].standardName = $('input#' + variables[v].id + '.standardName')[0].value;
		variables[v].longName = $('input#' + variables[v].id + '.longName')[0].value;
		variables[v].units = $('input#' + variables[v].id + '.units')[0].value;
		variables[v].bareGround = $('input#' + variables[v].id + '.bareGround')[0].checked;
		
		for (var i = 0; i < 5; i++)
			for (var j = 0; j < 4; j++) {
				var id = variables[v].id * 1000 + i * 10 + j;
				if ($('input#' + id + '.checkbox')[0].checked) {
					var freq = (Math.floor(id / 10)) % 10;
					var form = id % 10;
					var nic = $('input#' + id + '.nameInCode')[0].value;

					var variant = new Variant(variables[v].id, freq, form, nic);
					variants.push(variant);
				}
			}
	}
}

function generateVariants() {
	$('#outputViewer').empty();
	$('#outputViewer').append($('<p/>')
			.append($('<label for="filename"\>').width('100px').text('File Name:'))
			.append($('<input type="filename" id="filename">').attr('value', 'output.xml').change(
					function() {
						$('#downloadButton').attr('download', $('input#filename')[0].value);
					}))
			.append('&nbsp;&nbsp;')
			.append($('<a/>').attr('id', 'downloadButton').text('Download XML').attr('download', 'output.xml')).button());
	$('#outputViewer')
	.prepend(
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
			case 0: $freq.html('annually'); break;
			case 1: $freq.html('monthly'); break;
			case 2: $freq.html('daily'); break;
			case 3: $freq.html('6 hourly'); break;
			case 4: $freq.html('halfhourly'); break;
			}

			switch (form) {
			case 0: $form.html('grid'); break;
			case 1: $form.html('pft'); break;
			case 2: $form.html('tile'); break;
			case 3: $form.html('layer'); break;
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
