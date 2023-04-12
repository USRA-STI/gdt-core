.. _core-collection:
.. |DataCollection| replace:: :class:`~gdt.core.collection.DataCollection`

*********************************************
Data Collections (:mod:`gdt.core.collection`)
*********************************************

.. currentmodule:: gdt.core.collection

The DataCollection Class
========================
Often one has multiple data objects and would like to easily perform the same
actions over the collection of objects.  The |DataCollection| class is designed
specifically for this purpose.

Examples
--------
As an example, we will create a trivial class to demonstrate how the 
DataCollection works.

    >>> class Trivial():
    >>>     def __init__(self, filename, attribute):
    >>>         self._filename = filename
    >>>         self._attribute = attribute
    >>>
    >>>     @property
    >>>     def filename(self):
    >>>         return self._filename
    >>>
    >>>     def get_attribute(self):
    >>>         return self._attribute
    >>>
    >>>     def set_attribute(self, attr):
    >>>         self._attribute = attr

We can pretend our Trivial class contains data from some file, including some
attribute that we can retrieve and set.  If we have a bunch of these, we can use
a DataCollection:

    >>> from gdt.core.collection import DataCollection
    >>> triv1 = Trivial('data1.fit', 'yolo')
    >>> triv2 = Trivial('data2.fit', 'lol')
    >>> triv3 = Trivial('data3.fit', 'iykyk')
    >>> collection = DataCollection.from_list([triv1, triv2, triv3])
    >>> collection
    <DataCollection: 3 Trivial objects>
    
You can retrieve the number of objects in the collection using ``len()`` and
you can also retrieve the names associated with the objects in the collection:

    >>> len(collection)
    3
    >>> collection.items
    ['data1.fit', 'data2.fit', 'data3.fit']

The DataCollection automatically associates a ``filename`` attribute with the
object, but different names can be assigned when creating a the DataCollection
object.  For example:

    >>> collection = DataCollection.from_list([triv1, triv2, triv3], 
    >>>                                        names=['name1', 'name2', 'name3'])
    >>> collection.items
    ['name1', 'name2', 'name3']

You can add new objects, as long as they are of the same type already in the
collection:

    >>> # cannot do this
    >>> collection.include('this is not a Trivial Object', 'whoops')
    TypeError: A DataCollection must contain like objects
    
    >>> # can do this
    >>> triv4 = Trivial('data4.fit', 'brb')
    >>> collection.include(triv4, name='name4')
    >>> collection.items
    ['name1', 'name2', 'name3', 'name4']

And you can remove items by their name:
    
    >>> collection.remove('name4')
    >>> collection.items
    ['name1', 'name2', 'name3']
    
Of course, the most powerful aspect of the DataCollection is that it exposes
all attributes and methods of the objects in the collection so that you
only need to call them once.  For example, if we want to retrieve the filenames
for all of our objects:

    >>> collection.filename()
    ['data1.fit', 'data2.fit', 'data3.fit']

Note that even though ``filename`` is an attribute of the Trivial class, it is
called as a method here.  Similarly, if we want to call a function to retrieve
information from the collected objects:

    >>> collection.get_attribute()
    ['yolo', 'lol', 'iykyk']
    
You can similarly call a function that takes arguments and keywords or even
sets properties.  One thing to remember is that all arguments and keywords are
applied to all objects in the collection:
    
    >>> collection.set_attribute('same')
    >>> collection.get_attribute()
    ['same', 'same', 'same']
    
If you want to use different arguments for each object, you can iterate:
    
    >>> attrs = ['veni', 'vidi', 'vici']
    >>> for item, attr in zip(collection, attrs):
    >>>     item.set_attribute(attr)
    >>> collection.get_attribute()
    ['veni', 'vidi', 'vici']

Finally, you can retrieve an object from the collection by its name or return
the contents of the collection as a list:

    >>> collection.get_item('name2')
    <__main__.Trivial at 0x106100e90
    
    >>> collection.to_list()
    [<__main__.Trivial at 0x1060ef1d0>,
     <__main__.Trivial at 0x106100e90>,
     <__main__.Trivial at 0x106100710>]



Reference/API
=============

.. automodapi:: gdt.core.collection
   :inherited-members:


