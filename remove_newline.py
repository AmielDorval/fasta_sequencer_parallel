def remove_newlines(input_file, output_file):
    with open(input_file, 'r') as infile:
        with open(output_file, 'w') as outfile:
            for line in infile:
                stripped_line = line.rstrip('\n')
                outfile.write(stripped_line)

input_file = 'expected_5k.txt'  # Replace with the path to your input file
output_file = 'expected_5k_noNewline.txt'  # Replace with the path to your output file
remove_newlines(input_file, output_file)
