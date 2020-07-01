

class ExampleClass():
    def __init__(self, val):
        print('creating ExampleClass')
        self.value = val

    def hello_in_class(self):
        print('hello from ExampleClass')
        
    def return_value_plus(self, add):
        return self.value + add
        
