# Define the file paths
base_path = 'simulationSetups/deeponet/'
input_file_path = base_path + 'tdesign_60_points.txt'
output_file_path = base_path + 'setup_test_srcpos_fixed/cube6x6x6_src27/receivers_radius1m_origo333_sphere60.txt'

# Function to read data from the input file, add a constant, and write to the output file
def process_data(input_file, output_file, constant):
    with open(input_file, 'r') as f:
        # Read lines from the input file and process each line
        lines = f.readlines()
        processed_data = []

        for line in lines:
            # Split each line into a list of floating-point numbers
            values = [float(x) for x in line.split()]
            
            # Add the constant to each value
            updated_values = [x + constant for x in values]

            # Append the updated values to the processed data list
            processed_data.append(updated_values)

    with open(output_file, 'w') as f:
        # Write the processed data to the output file
        for values in processed_data:
            formatted_values = ['{:.4f}'.format(x) for x in values]
            # Convert the list of values to a space-separated string
            line = ' '.join(formatted_values)
            f.write(line + '\n')

constant_to_add = 3.0
process_data(input_file_path, output_file_path, constant_to_add)

print(f"Data has been processed and written to {output_file_path}.")