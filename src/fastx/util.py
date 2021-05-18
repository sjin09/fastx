from collections import defaultdict

class NestedDefaultDict(defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

    def __repr__(self):
        return repr(dict(self))


def chunkstring(string):
    chunks = [
        string[i : i + 60] for i in range(0, len(string), 60)
    ]
    return chunks

