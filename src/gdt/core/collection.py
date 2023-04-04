# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
# Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland and Adam Goldstein
#               Universities Space Research Association
#               Science and Technology Institute
#               https://sti.usra.edu
#
# Developed by: Daniel Kocevski
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#
import re
from functools import partial
from collections import OrderedDict

__all__ = ['DataCollection']

class DataCollection():
    """A container for a collection of like data objects, such as a collection
    of :class:`~.gdt.core.phaii.Phaii` objects.  This class exposes the individual objects'
    attributes, and it exposes the methods of the object class so that methods 
    can be called on the collection as a whole.  For that reason, each object in 
    the DataCollection must be of the same type, otherwise an error is raised.  
    The type of the collection is set by the first object inserted into the 
    collection and the collection type is immutable thereafter.
    
    Objects are stored in the collection in the order they are added.
    
    The number of items in the collection can be retrieved by ``len()`` and
    one can iterate over the items:: 
        
        [data_item for data_item in DataCollection]
    
    In addition to the DataCollection methods, all of the individual object 
    attributes and methods are exposed, and they become methods of the 
    DataCollection.  Note that individual object attributes become *methods* 
    i.e. if you have an item attribute called item.name, then the corresponding 
    DataCollection method would be item.name().
    """
    def __init__(self):
        self._data_dict = OrderedDict()
        self._type = None

    def __iter__(self):
        for item in self._data_dict.values():
            yield item

    def __len__(self):
        return len(self._data_dict)

    @property
    def items(self):
        """(list): The names of the items in the DataCollection"""
        return list(self._data_dict.keys())

    @property
    def types(self):
        """(str): The type of the objects in the DataCollection"""
        return self._type

    @classmethod
    def from_list(cls, data_list, names=None):
        """Given a list of objects and optionally a list of corresponding names, 
        create a new DataCollection. 
        
        Args:
            data_list (list of :obj:`objects`): 
                The list of objects to be in the collection
            names (list of :obj:`str`, optional):  
                The list of corresponding names to the objects.  If not set, 
                will try to retrieve a name from object.filename (assuming it's 
                a data object). If that fails, each item will be named 
                ambiguously 'item1', 'item2', etc.
        
        Returns                
            (:class:`~gdt.core.collection.DataCollection`)
        """
        obj = cls()

        # set the names
        if names is not None:
            if len(names) != len(data_list):
                raise ValueError('Names list must be same size as data list')
        else:
            names = [None] * len(data_list)

        # include the objects
        for data_item, name in zip(data_list, names):
            obj.include(data_item, name=name)

        return obj

    def get_item(self, item_name):
        """Retrieve an object from the DataCollection by name
        
        Args:
            item_name (str): The name of the item to retrieve
        
        Returns:
            (:obj:`object`)
        """
        return self._data_dict[item_name]

    def include(self, data_item, name=None):
        """Insert an object into the collection.  The first item inserted will 
        set the immutable type.
        
        Args:
            data_item (:obj:`object`): A data object to include
            name (str, optional): 
                An optional corresponding name.  If not set, will try to 
                retrieve a name from object.filename (assuming it's a data 
                object). If that fails, each item will be named ambiguously 
                'item1', 'item2', etc.
        """
        # if this is the first item inserted, set the type of the Collection
        # and expose the attributes and methods of the object
        if len(self) == 0:
            self._type = type(data_item)
            dir = [key for key in data_item.__dir__() if
                   not re.match('_.', key)]
            for key in dir:
                setattr(self, key, partial(self._method_call, key))   
        else:
            # otherwise, ensure that each object inserted is of the same type
            if type(data_item) != self._type:
                raise TypeError('A DataCollection must contain like objects')

        # insert with user-defined name
        if name is not None:
            self._data_dict[name] = data_item
        else:
            # or try to insert using filename attribute
            try:
                self._data_dict[str(data_item.filename)] = data_item
            # otherwise default to ambiguity
            except AttributeError:
                self._data_dict['item{}'.format(len(self) + 1)] = data_item

    def remove(self, item_name):
        """Remove an object from the collection given the name 
        
        Args:
            item_name (str): The name of the item to remove
        """
        self._data_dict.pop(item_name)

    def to_list(self):
        """Return the objects contained in the DataCollection as a list.
        
        Returns:
            (list of :obj:`objects`)
        """
        return [self.get_item(name) for name in self.items]

    def _enforce_type(self, data_item):
        if not isinstance(data_item, self._type) and self._type is not None:
            raise TypeError(
                'Incorrect data item for {}'.format(self.__class__.__name__))

    def _method_call(self, method_name, *args, **kwargs):
        """This is the wrapper for the exposed attribute and method calls.  
        Applies method_name over all items in the DataCollection
        
        Args:
            method_name (str): The name of the method or attribute
            *args: Additional arguments to be passed to the method
            **kwargs: Additional keyword arguments to be passed to the method
        
        Returns:
            None or list: If not None, will return the results from all 
            objects in the list     
       """
        # get the attributes/methods for each item
        refs = [getattr(obj, method_name) for obj in self._data_dict.values()]

        # if method_name is a method, then it will be callable
        if callable(refs[0]):
            res = [getattr(obj, method_name)(*args, **kwargs)
                   for obj in self._data_dict.values()]
        # otherwise, method_name will not be callable if it is an attribute
        else:
            # we are setting an attribute    
            if len(args) != 0:
                res = [setattr(obj, method_name, *args)
                       for obj in self._data_dict.values()]
            # we are retrieving an attribute
            else:
                res = refs

        if res[0] is not None:
            return res

    def __repr__(self):
        if len(self) > 0:
            _type = self.to_list()[0].__class__.__name__
        else:
            _type = ''
        s = '<{0}: {1} {2} objects>'.format(self.__class__.__name__,
                                            len(self), _type)
        return s
