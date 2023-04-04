.. _core-detectors:
.. |Detectors| replace:: :class:`~gdt.core.detector.Detectors`

*************************************************
Instrument Detectors  (:mod:`gdt.core.detectors`)
*************************************************

Introduction
============
Sometimes a gamma-ray instrument is composed of many individual detectors, and
those detectors may have different viewing directions.  The |Detectors| class
provides a simple way to store the detector names, indexing, and pointings so
that they can be retrieved for a variety of uses.

For Developers:
===============
The |Detectors| class cannot be used directly; it must be sub-classed with 
the relevant information for each detector.  As an example, we will create a
new class called ``MyDetectors`` that inherits the Detectors class, and we
will add a few detector definitions:

  >>> from gdt.core.detector import Detectors
  >>> class MyDetectors(Detectors):
  >>>     det0 = ('Detector0', 0, 0.0, 15.0)
  >>>     det1 = ('Detector1', 1, 15.0, 30.0)
  >>>     det2 = ('Detector2', 2, 30.0, 45.0)
  >>>     det3 = ('Detector3', 3, 45.0, 60.0)

In the example, we sub-classed Detectors and added four detector definitions,
assigned each one to a class variable.  The class variable name for each 
detector should be a short-hand name of the detector.  In each definition, there
are four values: the full name of the detector, the index of the detector, 
the azimuth angle of the detector normal in spacecraft coordinates, and the
zenith angle of the detector normal in spacecraft coordinates.

Examples
--------

Based on the use case, you may want to retrieve the detector information using
the short name, full name, or index number.  There are four distinct ways to
retrieve a detector from our Detectors class:

  >>> MyDetectors.det0
  <MyDetectors: det0>
  >>> MyDetectors.from_str('det1')
  <MyDetectors: det1>
  >>> MyDetectors.from_full_name('Detector2')
  <MyDetectors: det2>
  >>> MyDetectors.from_num(3)
  <MyDetectors: det3>
  
You can also iterate over the detectors in the Detectors class:

  >>> dets = [det for det in MyDetectors]
  >>> dets
  [<MyDetectors: det0>,
   <MyDetectors: det1>,
   <MyDetectors: det2>,
   <MyDetectors: det3>]

For a detector, you can retrieve any of the properties in the detector 
definition:
    
  >>> MyDetectors.det0.azimuth
  0.0
  >>> MyDetectors.det0.zenith
  15.0
  >>> MyDetectors.det1.pointing()
  (15.0, 30.0)
  >>> [det.full_name for det in MyDetectors]
  ['Detector0', 'Detector1', 'Detector2', 'Detector3']

Reference/API
=============

.. automodapi:: gdt.core.detector
   :inherited-members:
