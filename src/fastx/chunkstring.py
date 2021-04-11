def chunkstring(string, string_length=50):
    chunks = [
        string[i : i + string_length] for i in range(0, len(string), string_length)
    ]
    return chunks

