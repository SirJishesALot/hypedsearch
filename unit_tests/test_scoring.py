class test_utils(unittest.TestCase):
    def setUp(self):
        self.tosortkey = [
            {'a': 'b', 'b': 2}, 
            {'a': 'e', 'b': 18},
            {'a': 'd', 'b': 11},
            {'a': 'a', 'b': -1}, 
            {'a': 'c', 'b': 3}    
        ]
        self.tosortindex = [
            ('b', 2), ('e', 18), ('d', 11), ('a', -1), ('c', 3)
        ]