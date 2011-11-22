
class Links(dict):
	""" An extended dictionary class 
			that handles tuple keys nicely,
			by ordering the tuple calls """

	def __init__(self,info=None,**kwargs):
		super(Links,self).__init__()
		if info:
			if hasattr(info,'iteritems'):
				for key,value in info.iteritems():
					self.__setitem__(key,value)
			elif iterable(info):
				for key,value in info:
					self.__setitem__(key,value)
		for key,value in kwargs.iteritems():
				self.__setitem__(key,value)
	def order_first_tuple(func):
		def orderedversion(self,first,*args,**kwargs):
			assert hasattr(first,"__len__"), "Key must be iterable"
			assert len(first)==2, "Key must be 2-tuple"
			ordered_first = tuple(sorted(first))
			return func(self,ordered_first,*args,**kwargs)
		return orderedversion
	@order_first_tuple
	def __setitem__(self,key,val):
		return super(Links,self).__setitem__(key,val)

	def __getitem__(self,first,*args):
		if 
		return super(Links,self).__getitem__(*args)
	@order_first_tuple
	def __delitem__(self,*args):
		return super(Links,self).__delitem__(*args)
	@order_first_tuple
	def __contains__(self,*args):
		return super(Links,self).__contains__(*args)