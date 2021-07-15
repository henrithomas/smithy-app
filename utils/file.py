def read_large_file(file_object):
    while True:
        data = file_object.readline()
        if not data:
            break
        yield data