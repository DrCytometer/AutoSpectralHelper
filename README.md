# AutoSpectralHelper: Interactive Control File Builder

AutoSpectralHelper is a Shiny application for building and validating the spectral flow cytometry control
files for use in AutoSpectral. It provides an interactive interface for selecting controls, marking 
unstained samples, editing metadata, and exporting a finalized control file.


## Features

- Select cytometer configuration
- Load control files from a directory
- Interactive rHandsontable control editor
- Mark/unmark unstained controls for use as universal negatives
- Reorder controls (helpful for plotting fluorophore spectra later)
- Auto-fill functions
- Export validated control files


## Installation

Download the `app.r` file.

While you can manually install the dependencies, this occurs automatically when you first run the app.

You will need AutoSpectral version 0.8.7 or higher. Just re-install it if you aren't sure.
```r
devtools::install_github("DrCytometer/AutoSpectral")
```


## Using the shiny helper

Place a copy of the app in the folder where you want to work. The app will look for your single-stained 
control files in this folder and its contents, with your assistance.

To launch the app, click the green "Run App" button on the upper right of the app file.

For more details and instructions on creating the control file, visit the help page on AutoSpectral:
https://drcytometer.github.io/AutoSpectral/articles/Control_File_example.html 


## Requirements

R >= 4.0
Shiny >= 1.7
AutoSpectral >= 0.8.7


## Issues

If you encounter problems, please open an issue:
https://github.com/DrCytometer/AutoSpectralHelper/issues


## License

AGPL (>= 3)
