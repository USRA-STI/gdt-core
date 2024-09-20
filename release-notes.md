# Release Notes for Gamma-ray Data Tools: Core Library
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
- Fits files are now support context management and perform header verification.
