


class test:

    def __init__(self):

        self.a = 1
        self.b = 2
        self.func1()

    def func1(self):
        self.a = 10
        self.b = 20
        self.c = 3

    def func2(self):
        print(self.__dict__)


t = test()
t.func2()
