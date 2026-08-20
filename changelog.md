# Changelog

All notable changes to this project will be documented in this file.

## [0.3.1] - Adding BPMN Colour According to Skip Probabilities - #3

### Added

* Added a visualisation tool for calculated skip probabilities.
* Added functionality to colour BPMN activities according to their calculated skip probabilities.
* Added support for generating a visually annotated `.bpmn` file to make skip probabilities easier to interpret.

### Purpose

This feature provides a visual representation of skip probabilities directly within the BPMN model, allowing activities with different skip probabilities to be identified more easily. 

### Implementation
The colours added depend upon tools such as bpmn.io to view the model. 
