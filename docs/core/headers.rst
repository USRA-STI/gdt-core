.. _core-headers:
.. |Header| replace:: :class:`~gdt.core.headers.Header`
.. |FileHeaders| replace:: :class:`~gdt.core.headers.FileHeaders`

*******************************************
Data File Headers (:mod:`gdt.core.headers`)
*******************************************

Introduction
============
One of the most common standardized file formats in astronomy is the FITS
format.  In this format, files can have one or more extension, where each
extension contains a header and optionally some data.  The standard is that 
every file has a ``Primary`` header containing high-level metadata about the
instrument, date, and other properties pertaining to how the data in the file
was collected.  Other extensions containing data will also have headers with
metadata relevant to the data in that extension.  The GDT provides two classes
to assist in reading, modifying, or creating the metadata information in the
FITS file headers.

The Header Class
================
Astropy already provides an API for reading, modifying, and creating a FITS
header, but we want to impose some important limitations for quality control,
so the |Header| class in the GDT is a sub-class of the Astropy Header class.  
Specifically, the GDT serves as a way to create FITS headers that follow a 
prescribed template with values that are *strictly* typed and disallow 
addition of new cards.  These restrictions are important for maintaining high
quality of header information when creating and modifying file headers for 
production and public use.

For Developers:
---------------
To use the |Header| class, it should be sub-classed with two class variables:
  * ``name`` (str): The extension name for the header
  * ``keywords`` (list): A list of header cards, which are 3-tuples containing
                         the keyword name, default value, and comment string

The list of cards in the ``keywords`` list will exhaustively define the 
cards that will be in the header.  No additional cards can be added once the 
header is instantiated.  The default value for each card defines the data type
of the value. For example, a string value requires that if the value is updated,
it must be a string.

We will make a header definition called ``MyHeader`` by sub-classing |Header|:

  >>> from gdt.core.headers import Header
  >>> class MyHeader(Header):
  >>>     name = 'PRIMARY'
  >>>     keywords = [('STRING', 'hello', 'A string value'),
  >>>                 ('INT', 1, 'An integer value'),
  >>>                 ('FLOAT', 5.7, 'A float value'),
  >>>                 ('BOOL', True, 'A Boolean value')]
  
This is our definition for our new header, which is the primary header.  We 
have defined four cards, the first one with keyword 'STRING', that has a default
value of 'hello,' and a brief comment describing the value.

Examples
--------

Here is what it looks like when we instantiate ``MyHeader``:

  >>> header = MyHeader()
  >>> header
  STRING  = 'hello   '           / A string value                                 
  INT     =                    1 / An integer value                               
  FLOAT   =                  5.7 / A float value                                  
  BOOL    =                    T / A Boolean value                                

You can access the keywords as normal, or even update their values, but keep in
mind that the type for each keyword is fixed upon initialization and cannot be
changed:

  >>> header['STRING']
  'hello'
  >>> # update keyword value
  >>> header['INT'] = 22
  >>> header['INT']
  22
  >>> # update keyword value; type is enforced
  >>> header['INT'] = '50'
  >>> header['INT']
  50
  >>> # attempted update keyword with invalid value
  >>> header['INT'] = 'oops'
  TypeError: Value for INT is of incorrect type

Finally, you can set values for the keywords upon initialization of the header
by defining the associated kwargs:

  >>> header = MyHeader(string='goodbye', int=250, float=42.7, bool=False)
  >>> header
  STRING  = 'goodbye '           / A string value                                 
  INT     =                  250 / An integer value                               
  FLOAT   =                 42.7 / A float value                                  
  BOOL    =                    F / A Boolean value                                

Note that if you have a keyword in the header that has a hyphen 
(e.g. 'DATE-OBS'), the associated kwarg will replace the hyphen with an 
underscore (e.g. date_obs).


The FileHeaders Class
=====================
Typically a FITS file with have multiple extensions, and therefore multiple 
headers.  The |FileHeaders| class helps manage multiple headers belonging to
a single file, serving multiple purposes: the ability to efficiently and simply
create a default set of headers for a file, manage the update of values across
multiple headers in a file, and acting as a type of verification for headers
read from a file.

For Developers:
---------------
The FileHeaders class should not be instantiated directly but instead should
be sub-classed with one required class variable:

  * ``_header_templates`` (list): A list of |Header| objects that belong to the file.

This class variable is used to build the expected default set of headers and
use the header definitions as a verification against headers that are read in
from a file.

Let us use our primary header from the example in the previous section and add
another header:

  >>> class MySecondHeader(Header):
  >>>     name = 'SECONDARY'
  >>>     keywords = [('ONE_KEY', 0, 'A keyword'),
  >>>                 ('DATE', '', 'The date'),
  >>>                 ('COMMENT', 'foobar', 'A comment')]

And then we can define our |FileHeaders| class called ``MyFileHeaders``:

  >>> from gdt.core.headers import FileHeaders
  >>> class MyFileHeaders(FileHeaders):
  >>>     _header_templates = [MyHeader(), MySecondHeader()]

That's it.  And now we create the object:

  >>> file_headers = MyFileHeaders()
  >>> file_headers
  <MyFileHeaders: 2 headers>

Examples
--------
You can retrieve the list of extension names and extract the headers by index
or extension name:

  >>> file_headers.keys()
  ['PRIMARY', 'SECONDARY']
  >>> file_headers['PRIMARY']
  STRING  = 'hello   '           / A string value                                 
  INT     =                    1 / An integer value                               
  FLOAT   =                  5.7 / A float value                                  
  BOOL    =                    T / A Boolean value
  >>> file_headers[1]
  ONE_KEY = '' / A keyword                                                        
  DATE    = '2022-04-23T15:11:28' / The date                                      
  COMMENT foobar                                                                  

You will notice the ``DATE`` keyword in the second header is filled in even 
though we specified a null string as default.  The ``DATE`` keyword is a special
keyword that specifies the creation date of the header.  If the header is part
of a file to be written to disk, you should call the ``update()`` method to 
update the value to the current date and time:

  >>> file_headers.update()
  >>> file_headers[1]['DATE']
  '2022-04-23T17:11:05'

The FileHeaders class can act as a header verifier for reading in files. If 
there are missing keywords compared to what is expected in the header template,
or if the values are not the type expected, an exception will be raised. If 
there are extraneous keywords in the header compared to what is expected, no
exception will be raised, but the extraneous keywords will be ignored.

To test this, let's make a similar header to ``MyHeader``, but with the ``BOOL``
keyword missing:

  >>> class MyPartialHeader(Header):
  >>>     name = 'PRIMARY'
  >>>     keywords = [('STRING', 'hello', 'A string value'),
  >>>                 ('INT', 1, 'An integer value'),
  >>>                 ('FLOAT', 5.7, 'A float value')]

And now we can read in a list of headers to ``MyFileHeaders``. It will verify
the headers are consistent with the expected header definitions, and if so,
will populate the template headers with the header information:

  >>> # the list of headers can be from a file read from disk
  >>> file_headers = MyFileHeaders.from_headers([MyPartialHeader(), MySecondHeader()])
  KeyError: "Keyword 'BOOL' not found."
  >>> file_headers = MyFileHeaders.from_headers([MyHeader(), MySecondHeader()])
  >>> file_headers
  <MyFileHeaders: 2 headers>
  

Reference/API
=============

.. automodapi:: gdt.core.headers
   :inherited-members:
