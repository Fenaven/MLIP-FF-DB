import hashlib

def print_hash(input_string):
    hash_object = hashlib.sha256()
    
    hash_object.update(input_string.encode('utf-8'))
    
    return hash_object.hexdigest()

# Example usage:
input_string = "orca 1"
print_hash(input_string)
print(type(print_hash(input_string)))