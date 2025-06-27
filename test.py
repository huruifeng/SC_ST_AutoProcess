import random
import string
import os
import time


# Define file paths
file1_path = "file1.txt"
file2_path = "file2.txt"
file3_path = "file3.txt"
file4_path = "file4.txt"

# Helper function to write a file
def generate_file(path, col1_generator, col2_range, num_lines):
    with open(path, "w") as f:
        for _ in range(num_lines):
            col1 = next(col1_generator)
            col2 = random.randint(*col2_range)
            f.write(f"{col1}\t{col2}\n")

def generate_file_4cols(path, col1_generator, col2_range, num_lines):
    with open(path, "w") as f:
        for _ in range(num_lines):
            col1 = next(col1_generator)
            col2 = random.randint(*col2_range)
            col3 = random.randint(*col2_range)
            col4 = random.randint(*col2_range)
            f.write(f"{col1}\t{col2}\t{col3}\t{col4}\n")

# Generator for file1 column1: A000001 to A100000
def col1_gen_file1():
    while True:
        yield f"AA{random.randint(1, 100000):07d}"

# Generator for file2 column1: A01 to A99
def col1_gen_file2():
    while True:
        yield f"A{random.randint(1, 99):02d}"

    

# Generate files
generate_file(file1_path, col1_gen_file1(), (10000, 99999), 50_000_000)
generate_file(file3_path, col1_gen_file1(), (10000, 99999), 25_000_000)
generate_file(file2_path, col1_gen_file2(), (10, 99), 50_000_000)
generate_file_4cols(file4_path, col1_gen_file1(),  (10000, 99999), 50_000_000)

# Function to test read time
def test_read_time(file_path):
    start_time = time.time()
    with open(file_path, "r") as f:
        lines = f.readlines()
    end_time = time.time()
    return end_time - start_time, len(lines)

# Measure read time for both files
time_file1, lines_file1 = test_read_time(file1_path)
time_file2, lines_file2 = test_read_time(file2_path)
time_file3, lines_file3 = test_read_time(file3_path)
time_file4, lines_file4 = test_read_time(file4_path)

print("File1 Time:",time_file1, "Lines:",lines_file1)
print("File2 Time:",time_file2, "Lines:",lines_file2)
print("File3 Time:",time_file3, "Lines:",lines_file3)
print("File4 Time:",time_file4, "Lines:",lines_file4)

