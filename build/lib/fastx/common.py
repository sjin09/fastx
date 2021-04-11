from collections import defaultdict

class NestedDefaultDict(defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

    def __repr__(self):
        return repr(dict(self))


def chunkstring(string, string_length=50):
    chunks = [
        string[i : i + string_length] for i in range(0, len(string), string_length)
    ]
    return chunks

