# rcaiman project

This directory contains a project of the rcaiman suite.

Projects are managed with **rcaiman-projects**, a desktop application for 
organizing canopy image analysis workflows.

**Project homepage**

https://github.com/GastonMauroDiaz/rcaiman-projects

## Working with this project

The application automatically discovers images, point sets and other project
resources using documented folder and file naming conventions.

You are free to add your own folders and files to this project whenever needed.

Standard project folders and file names, however, have semantic meaning. Renaming or
reorganizing them may prevent the application from automatically discovering 
related resources.

## Design philosophy

An rcaiman project is primarily a directory structure.

The filesystem is the source of truth.

Rather than relying on databases or proprietary project formats, rcaiman-projects 
interprets project contents through simple, documented conventions whenever possible.

This makes projects easy to inspect, archive, version, share and understand using 
standard filesystem tools.
