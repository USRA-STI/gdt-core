# Release Notes for Gamma-ray Data Tools: Core Library
## Version 2.2.1 (Released Nov 19, 2025)

This release included the following updates from pull requests:  

- Update HEASARC Catalog Link [#91](https://github.com/USRA-STI/gdt-core/pull/91)
- Added support for retrieving HEASARC data from AWS servers [#92](https://github.com/USRA-STI/gdt-core/pull/92)
- Fix the finder behavior when using the AWS protocol to search directories instead of files [#93](https://github.com/USRA-STI/gdt-core/pull/93)
- Added links to GitHub repo and related pages (issues, pull requests) to documents [#95](https://github.com/USRA-STI/gdt-core/pull/95)
- Moved setting the x-axis limits to after the axes were configured [#98](https://github.com/USRA-STI/gdt-core/pull/98)
- Fixes an issue with rebinned channelize spectrum having gaps in the step plots. [#102](https://github.com/USRA-STI/gdt-core/pull/102)
- Fix tte deadtime [#103](https://github.com/USRA-STI/gdt-core/pull/103)
- Fixes issues with HEASARC change directory [#107](https://github.com/USRA-STI/gdt-core/pull/107)

## Version 2.2.0 (Released Apr 17, 2025)

This release included the following updates from pull requests:  

- Bug fix for DataCollection where names is None [#62](https://github.com/USRA-STI/gdt-core/issues/62)
- Added detector propagation to rebin and resample. [#65](https://github.com/USRA-STI/gdt-core/pull/65)
- Added functions to create FITS headers with floating-point values. [#69](https://github.com/USRA-STI/gdt-core/pull/69)
- Updated calls to matplotlib and fixed some headers. [#71](https://github.com/USRA-STI/gdt-core/pull/71)
- Renamed exponential_card to scientific_card and added some options. [#72](https://github.com/USRA-STI/gdt-core/pull/72)
- Reduced polygon artifacts for maps with fragmented contours. [#77](https://github.com/USRA-STI/gdt-core/pull/77)
- Fixed creating and plotting effective area on the sky. [#81](https://github.com/USRA-STI/gdt-core/pull/81)
- Fixed 'Pha.energy_range' when creating from Phaii.to_pha(). [#82](https://github.com/USRA-STI/gdt-core/pull/82)
- Improved the creation of PHA objects with 'valid_channels'. [#83](https://github.com/USRA-STI/gdt-core/pull/83)
- Removed the need for using a FILENAME in the header of a PHA file. [#84](https://github.com/USRA-STI/gdt-core/pull/84)
- Added translation function to Range. [#85](https://github.com/USRA-STI/gdt-core/pull/85)
- Fixed more calls to matplotlib. [#86](https://github.com/USRA-STI/gdt-core/pull/86)
- Added support for trigger detection. [#90](https://github.com/USRA-STI/gdt-core/pull/90)

## Version 2.1.0 (Released Aug 16, 2024)  

This release included the following updates from pull requests:  

- Bug fix in creating annulus [#58](https://github.com/USRA-STI/gdt-core/pull/58)
- Re-arrangement of HEASARC HTTPS/FTP support [#57](https://github.com/USRA-STI/gdt-core/pull/57)
- Updated calls to Matplotlib [#56](https://github.com/USRA-STI/gdt-core/pull/56)
- Remove Need to Run "gdt-data init" manually [#52](https://github.com/USRA-STI/gdt-core/pull/52)
- SpectralFitter asymmetric errors fix [#55](https://github.com/USRA-STI/gdt-core/pull/55)
- Modified 'icrs_to_spacecraft' to remove 'Nan' vale of 'el'. [#45](https://github.com/USRA-STI/gdt-core/pull/45) 
- avoid list conversion by creating an empty numpy array to get a >2x speedup when opening TTE files [#44](https://github.com/USRA-STI/gdt-core/pull/44)
- Documentation on Contributions [#42](https://github.com/USRA-STI/gdt-core/pull/42)
- Duration cleanup [#40](https://github.com/USRA-STI/gdt-core/pull/40)
- Skyplot heatmap fix [#38](https://github.com/USRA-STI/gdt-core/pull/38)
- Fix for SuperFunction `names` attribute [#35](https://github.com/USRA-STI/gdt-core/pull/35)
- A more generalized fix for the attribute error fix for _repr_html(). [#37](https://github.com/USRA-STI/gdt-core/pull/37)
- Make DataCollection more pythonic with [] access [#11](https://github.com/USRA-STI/gdt-core/pull/11)
- Ensure FitsFileContextManager._repr_html_ runs when used with GbmHealPix.open [#33](https://github.com/USRA-STI/gdt-core/pull/33)
- Remove try/except [#32](https://github.com/USRA-STI/gdt-core/pull/32)
- Improve PG-stat stability [#29](https://github.com/USRA-STI/gdt-core/pull/29)
- Headers warning [#25](https://github.com/USRA-STI/gdt-core/pull/25)
- Pstat now handles zero-count bins [#27](https://github.com/USRA-STI/gdt-core/pull/27)  

## Version 2.0.4 (Released Feb 23, 2024)  
Note: Version was skipped due to a version number collision with PyPi.

This release included the following updates from pull requests:  

- Support for non-energy-calibrated data [#21](https://github.com/USRA-STI/gdt-core/pull/21)
- Added link to docs in README [#15](https://github.com/USRA-STI/gdt-core/pull/15)

## Version 2.0.2 (Released May 10, 2023)

This released fixed a typo in a type declaration that was causing incompatibility with Python 3.8


## Version 2.0.1 (Released Apr 13, 2023)

This released fixed minor issues that were discovered after release of 2.0.0

## Version 2.0.0 (Released Apr 12, 2023)

This is the initial release.

We started the version at 2.0.0 to indicate that this is our major API update from GBM Data Tools which is version 1.1.1

Changes from GBM Data Tools include:

- Making the library more generalized so that it can be used for other gamma-ray observatories besides GBM.
- Functions that is common to all gamma-ray observatories is released as Gamma-ray Data Tools: Core
- Functions that is specific to Fermi-GBM is released as Gamma-ray Data Tools: Fermi
- Mission Elapsed Time is stored within AstroPy Time which handles the conversion to other time systems.
- Spacecraft coordinates is now using Astropy Coordinate Frames.
- Quaternions use Numpy to vectorize operations and SciPy Rotation to translate between Quaternions and Direction Cosine Matrixes.
- Fits files now support context management and perform header verification.
