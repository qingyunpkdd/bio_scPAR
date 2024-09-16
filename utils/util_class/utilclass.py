# 在计算中存在大量的计算量，但是不一定每次都会用到，这时候就可以使用装饰器来实现懒加载

class lazyproperty:
    def __init__(self, func):
        self.func = func

    def __get__(self, instance, cls):
        if instance is None:
            return self
        else:
            value = self.func(instance)
            setattr(instance, self.func.__name__, value)
            return value