Output Variable Editor (OVE)
========

The Output Variable Editor (OVE) is a web interface that allows a user to configure the output variables for the CLASS-CTEM model. Using this interface the user can load an existing XML file and add/remove variables and variants. Once the changes are complete, the user may download an updated version of the XML configuration file.

# Opening Instructions
The *Output Variable Editor* folder is located in the *tools* folder of your local CLASS directory.
Open index.html in your favorite browser and edit the output variable descriptors using the web interface.

# Work flow
The various process steps are seen on the top of the page as tabs. The user first loads the existing XML document as a template, then optionally add/removes variables, then edits the variants associated with the variables, and finally generates the output XML that included the new variables and variant changes.
## Load the Existing XML
Browse to an existing XML file to add/remove variables and variants. Once the file has been selected, press the "Load XML" button to load the data.
## Add/Remove Variable
Fill in the **Standard Name**, **Short Name**, **Long Name**, **Units** fields and select a **Group** to classify the new variable. Select the **Include bare ground** check-box, if the variable includes bare ground. One may also use **Add variable group** to add a group or **Remove variable** to delete an existing variable.
## Variants editor
In the variants editor, the user can find a list of all of the existing variables. Click on one to find a detailed description of the variable as well as a list of potential variants. Check the box for particular variant and enter the desired **name in code**.
## Output XML
Once all the variable and variant changes are complete, proceed to the **Output XML** tab, where one can preview the newly generated XML document. Also, the website will trigger a download command, so the generated XML document will be pushed to the user for download.

# The XML Document Structure
The XML document follows a hierarchical structure which has the **variableSet** as a root node. This node has a number of **group** child nodes, as well as attributes such as *type, version and created*. The **group** nodes have a number of **variable** child nodes.

The **variable** nodes feature several attributes (e.g. *id, includeBareGround*) and fields (e.g. *standardName, longName, shortName,* etc. They also have a collection of **variant** child nodes. These variants are the individual instances of the variable, which can be identified by the unique **nameInCode** identifier. They also have a specific **timeFrequency** and **outputForm**.

# Validation / Schema
A relaxNG schema is provided for validation purposes.
In order to validate the XML Document, navigate to the schema folder and use the following command:

`java -jar jing.jar schema_1.0.xml path/myFile.xml`

where schema_1.0.xml is the current version of the schema and path/myFile.xml is the XML Document file of your choosing.

The schema will produce warnings and errors if there is any problem with the XML file, but it will not produce any output at all if everything is alright.

# Manual editing
There are certain functionalities that the present interface is lacking, like editing existing variable names, removing groups, editing groups and so on.
These tasks can be easily achieved by manually editing the XML input file.

**Do not manually edit the id fields of the variables, as that can cause serious problems with the variables**