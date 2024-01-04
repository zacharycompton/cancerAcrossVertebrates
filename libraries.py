import re
import os

def extract_libraries_from_r_script(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    library_lines = [line.strip() for line in lines if re.search(r'\blibrary\(', line)]

    return library_lines

def process_r_scripts_in_directory(directory_path):
    r_scripts = [file for file in os.listdir(directory_path) if file.endswith('.R')]

    all_libraries = []

    for script in r_scripts:
        script_path = os.path.join(directory_path, script)
        libraries = extract_libraries_from_r_script(script_path)
        all_libraries.extend(libraries)

    return all_libraries

if __name__ == "__main__":
    directory_path = "/Users/walkermellon/cancerAcrossVertebrates"
    libraries = process_r_scripts_in_directory(directory_path)

    print("Libraries found in R scripts:")
    for lib in libraries:
        print(lib)
