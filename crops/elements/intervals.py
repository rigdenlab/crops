from crops.about import __prog__, __description__, __author__, __date__, __version__

import copy

def _intervalise(subject):
    """Turns any integer or list of two integers into a :class:`~crops.elements.intervals.intinterval`.

    :param subject: Input parameter.
    :type subject: int, list, :class:`~crops.elements.intervals.intinterval`
    :raises TypeError: Argument must be an integer interval, an integer or a list of two integers.
    :return: A :class:`~crops.elements.intervals.intinterval` defined by the input.
    :rtype: :class:`~crops.elements.intervals.intinterval`

    """

    if isinstance(subject, intinterval):
        return subject

    raisemsg='Argument must be an intinterval, an integer or a list of two integers.'

    if isinstance(subject, list):
        if len(subject) != 2 and len(subject) != 0:
            raise TypeError(raisemsg)
        elif len(subject) == 2:
            for s in subject:
                if isinstance(s,float):
                    if s==int(s):
                        s=int(s)
                    else:
                        raise TypeError(raisemsg)
                if not isinstance(s,int):
                    raise TypeError(raisemsg)
            if subject[0]>subject[1]:
                subject=[[subject[1],subject[0]]]
            else:
                subject=[subject]

    elif isinstance(subject, int):
        subject=[[subject,subject]]
    else:
        raise TypeError(raisemsg)

    intervalised=intinterval(description='other',subint=subject)

    return intervalised

def _intervalize(subject):
    return _intervalise(subject)

class intinterval:
    """An integer interval object.
    The :class:`~crops.elements.intervals.intinterval` class represents a data structure to hold all
    non-connected sub-intervals making it up, extra information tags, and
    set operations for the intervals.
    It contains functions to organise intervals and make operations on them.

    :param description: An interval ID, defaults to 'intinterval'.
    :type description: str, optional
    :param subint: A list of two-integer lists, defaults to [].
    :type subint: list [int], optional
    :ivar tags: Useful information of the :class:`~crops.elements.intervals.intinterval`, including the default 'description'.
    :vartype tags: dict [any]
    :ivar subint: The list of sub-intervals in :class:`~crops.elements.intervals.intinterval`.
    :vartype subint: list [ list [int] ]

    :example:
    >>> from crops.elements import intervals as cint
    >>> my_interval = cint.intinterval('an interval')
    >>> my_interval2 = cint.intinterval('another interval')
    >>> my_interval = my_interval.union(15)
    >>> my_interval2 = my_interval2.union([2,7])
    >>> my_interval = my_interval.union(my_interval2)
    >>> print(my_interval)
    Integer interval object: (id="an interval", subsets="[[2,7],[15,15]]")
    >>> print(my_interval.intersection([3,19])) # Just printing, not setting the values
    Integer interval object: (id="an interval", subsets="[[3,7],[15,15]]")
    >>> my_interval.n_elements()
    7
    >>> my_interval.terminals()
    [2,15]
    >>> my_interval.contains(16)
    False
    >>> my_interval.contains(my_interval2)
    True
    >>> my_interval.contains([5,7])
    True

    """
    kind = 'Integer interval'
    __slots__= ['tags','subint']
    def __init__(self, description='intinterval', subint=None):
        self.tags={}
        self.tags['description']=description
        if subint is None:
            self.subint=[]
        else:
            self.subint = subint

    def __repr__(self):
        return self.kind+' object: (id="'+ self.description()+'", subsets='+str(self.subint)+')'
    def __iter__(self):
        return iter(self.subint)

    def description(self, newdescription=None):
        """Returns or modifies the 'description' tag.

        :param newdescription: If given, the value of 'description' in tags is replaced by 'newdescription', defaults to None
        :type newdescription: str, optional
        :return: If newdescription not given, the value of 'description' in tags is returned.
        :rtype: str

        """
        if 'description' not in self.tags:
            self.tags['description']=None
        if newdescription is None:
            return self.tags['description']
        else:
            self.tags['description']=newdescription

    def addtag(self, tag, value):
        """Adds a new tag to :class:`~crops.elements.intervals.intinterval`.

        :param tag: Key argument.
        :type tag: str
        :param value: Value argument.
        :type value: any
        :raises TypeError: If tag is not a string.

        """
        if not isinstance(tag,str):
            raise TypeError('Keys must be strings')
        self.tags[tag]=value

    def deltag(self, tag):
        """Deletes a tag from :class:`~crops.elements.intervals.intinterval`.

        :param tag: Key argument.
        :type tag: str
        :raises TypeError: If argument is not a string.

        """
        if not isinstance(tag,str):
            raise TypeError('Tags must be strings')
        if tag=='description':
            raise ValueError('Key "description" cannot be removed.')
        self.tags.remove(tag)

    def copy(self):

        return copy.copy(self)

    def deepcopy(self):

        return copy.deepcopy(self)

    def terminals(self, other=None):
        """Returns the first and last element in the :class:`~crops.elements.intervals.intinterval`.

        :param other: Defaults to None.
        :type other: int, list [int], :class:`~crops.elements.intervals.intinterval`, optional
        :return: A list of two integers indicating lower and higher limits of the interval (self if other is None, other otherwise).
            If input is an empty interval, this function will return an empty list.
        :rtype: list [int]

        """
        interval=self if other is None else _intervalise(other)
        if len(interval.subint)>0:
            return [interval.subint[0][0],interval.subint[-1][-1]]
        else:
            return []

    def n_elements(self,other=None):
        """ Returns the number of elements in the interval, subinterval, or any other set.

        :param other: Input argument, defaults to None.
        :type other: int, list [int], :class:`~crops.elements.intervals.intinterval`, optional
        :return: Number of elements in the interval (self if other is None, other otherwise).
        :rtype: int

        """

        n=0
        interval=self if other is None else _intervalise(other)
        for A in interval.subint:
            n += A[1]-A[0]+1

        return n

    def contains(self, other):
        """Checks if input interval is fully contained in self. B ⊂ A : A = self, B = other.

        :param other: Another interval.
        :type other: int, list [int], :class:`~crops.elements.intervals.intinterval`
        :return: Whether self contains other or not.
        :rtype: bool
        """
        # A ⊂ B : A = self, B = other
        other=_intervalise(other)
        if other.subint==[]:
            return True
        else:
            if self.subint==[]:
                return False
            if self.subint==other.subint:
                return True

        for B in other.subint:
            outcome=False
            for A in self.subint:
                if B[0]>=A[0] and B[1]<=A[1]:
                    outcome=True
                    break
            if not outcome:
                return False

        return True

    def union(self, other,newdesc=None):
        """A ⋃ B : A = self, B = other.

        :param other: Another interval.
        :type other: int, list [int], :class:`~crops.elements.intervals.intinterval`
        :param newdesc: If given, description of returned interval, defaults to None.
        :type newdesc: str, optional
        :return: The Union of both intervals.
        :rtype: :class:`~crops.elements.intervals.intinterval`

        """
        if newdesc is None:
            newdesc=self.description()
        # A ⋃ B : A = self, B = other
        other=_intervalise(other)
        if self.subint==[] or self.subint==other.subint:
            newint=self.deepcopy()
            newint.subint=copy.deepcopy(other.subint)
            return newint
        elif other.subint==[] and self.subint!=[]:
            newint=self.deepcopy()
            newint.description(newdesc)
            return newint

        newint=self.deepcopy()
        # Use of distributive property:
        # A = (A1 ⋃ A2); B = (B1 ⋃ B2)
        # A ⋃ B = (A ⋃ B1) ⋃ (A ⋃ B2)
        for B in other.subint:
            ingap=True
            newlist=[]
            for A in newint.subint:
                # other's left-hand side
                if B[0]>=A[0] and B[0]<=A[1]+1:
                    B[0]=A[0]
                    ingap=False
                # other's right-hand side
                if B[1]>=A[0]-1 and B[1]<=A[1]:
                    B[1]=A[1]
                    ingap=False
                if B[1]>A[1] and B[0]<A[0]:
                    ingap=False

            if ingap:
                for s in range(len(newint.subint)):
                    if s==0:
                        if B[1]<newint.subint[s][0]:
                            newlist.append(B)
                        newlist.append(newint.subint[s])
                    else:
                        if B[0]>newint.subint[s-1][1] and B[1]<newint.subint[s][0]:
                            newlist.append(B)
                            newlist.append(newint.subint[s])
                        else:
                            newlist.append(newint.subint[s])
                    if s==len(newint.subint)-1 and B[0]>newint.subint[s][1]:
                        newlist.append(B)
            else:
                for A in newint.subint:
                    if A[0]>=B[0] and A[1]<=B[1]:
                        if len(newlist)==0:
                            newlist.append(B)
                        else:
                            if newlist[-1]!=B:
                                newlist.append(B)
                    else:
                        newlist.append(A)

            newint.subint=copy.deepcopy(newlist)

        newint.description(newdesc)

        return newint

    def intersection(self,other,newdesc=None):
        """A ⋂ B : A = self, B = other.

        :param other: Another interval.
        :type other: int, list [int], :class:`~crops.elements.intervals.intinterval`
        :param newdesc: If given, description of returned interval, defaults to None.
        :type newdesc: str, optional
        :return: The Intersection of both intervals.
        :rtype: :class:`~crops.elements.intervals.intinterval`

        """
        if newdesc is None:
            newdesc=self.description()
        # A ⋂ B : A = self, B = other
        other=_intervalise(other)
        if self.subint==[] or self.subint==other.subint:
            newint=self.deepcopy()
            newint.description(newdesc)
            return newint
        elif other.subint==[] and self.subint!=[]:
            newint=self.deepcopy()
            newint.subint=[]
            newint.description(newdesc)
            return newint

        # Use of distributive property:
        # A = (A1 ⋃ A2); B = (B1 ⋃ B2)
        # A ⋂ B = (A1 ⋂ B1) ⋃ (A1 ⋂ B2) ⋃ (A2 ⋂ B1) ⋃ (A2 ⋂ B2)
        newint=self.deepcopy()
        newint.subint=[]
        for A in self.subint:
            for B in other.subint:
                partial_int=[]
                if B[0]>A[1] or B[1]<A[0]:
                    pass
                elif B[0]<A[0] and B[1]>=A[0] and B[1]<=A[1]:
                    partial_int=[A[0],B[1]]
                elif B[0]<A[0] and B[1]>A[1]:
                    partial_int=[A[0],A[1]]
                elif B[0]>=A[0] and B[1]<=A[1]:
                    partial_int=[B[0],B[1]]
                elif B[0]>=A[0] and B[0]<=A[1] and B[1]>A[1]:
                    partial_int=[B[0],A[1]]
                newint=newint.union(partial_int)

        newint.description(newdesc)

        return newint

    def subtract(self, other,newdesc=None):
        r"""A \\ B : A = self, B = other

        :param other: Another interval.
        :type other: int, list [int], :class:`~crops.elements.intervals.intinterval`
        :param newdesc: If given, description of returned interval, defaults to None.
        :type newdesc: str, optional
        :return: The result of subtracting other from self.
        :rtype: :class:`~crops.elements.intervals.intinterval`

        """
        if newdesc is None:
            newdesc=self.description()
        # A \ B : A = self, B = other
        other=_intervalise(other)
        if other.subint == []:
            newint=self.deepcopy()
            newint.description(newdesc)
            return newint
        elif other.subint == self.subint:
            newint=self.deepcopy()
            newint.subint=[]
            newint.description(newdesc)
            return newint

        newint=self.deepcopy()
        # Use of the following set algebra property:
        # A = (A1 ⋃ A2); B = (B1 ⋃ B2)
        # A \ B = [(A1 \ B1) \ B2] ⋃ [(A2 \ B1) \ B2]
        for A in newint.subint:
            for B in other.subint:
                if A!=[]:
                    if B[0]>A[1] or B[1]<A[0]:
                        pass
                    elif B[0]<A[0] and B[1]>=A[0] and B[1]<=A[1]:
                        if B[1]==A[1]:
                            A.clear()
                            break
                        else:
                            A[0]=B[1]+1
                    elif B[0]<=A[0] and B[1]>=A[1]:
                        A.clear()
                        break
                    elif B[0]>=A[0] and B[1]<=A[1] and B!=A:
                        if B[0]==A[0]:
                            A[0]=B[1]+1
                        elif B[1]==A[1]:
                            A[1]=B[0]-1
                        else:
                            newint.subint.insert(newint.subint.index(A)+1,[B[1]+1,A[1]])
                            A[1]=B[0]-1
                    elif B[0]>=A[0] and B[0]<=A[1] and B[1]>A[1]:
                        if B[0]==A[0]:
                            A.clear()
                            break
                        else:
                            A[1]=B[0]-1

        for i in reversed(range(len(newint.subint))):
            if newint.subint[i]==[]:
                newint.subint.pop(i)

        newint.description(newdesc)

        return newint

    def symdiff(self,other,newdesc=None):
        """A Δ B : A = self, B = other

        :param other: Another interval.
        :type other: int, list [int], :class:`~crops.elements.intervals.intinterval`
        :param newdesc: If given, description of returned interval, defaults to None.
        :type newdesc: str, optional
        :return: The Symmetric difference of both intervals.
        :rtype: :class:`~crops.elements.intervals.intinterval`

        """
        if newdesc is None:
            newdesc=self.description()
        # A Δ B : A = self, B = other
        other=_intervalise(other)
        if other.subint == []:
            newint=self.deepcopy()
            newint.description(newdesc)
            return newint
        elif self.subint==[] and other.subint!=self.subint:
            newint=self.deepcopy()
            newint.subint=copy.deepcopy(other.subint)
            newint.description(newdesc)
            return newint
        elif other.subint == self.subint:
            newint=self.deepcopy()
            newint.subint=[]
            newint.description(newdesc)
            return newint

        # Use of the following set algebra property:
        # A Δ B = (A \ B) ⋃ (B \ A)
        newint=self.subtract(other)
        newint=newint.union(other.subtract(self))

        newint.description(newdesc)

        return newint
