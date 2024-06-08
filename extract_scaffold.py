def description_process(input_file, output_file):
    with open(input_file, 'r') as input_file, open(output_file, 'w') as output_file:
        for line in input_file:
            if line.startswith('>'):
                parts = line.strip().split(',')
                output_file.write(f">{parts[0][-13:]}\n")
            else:
                output_file.write(line)

def genemark_process(input_file, output_file):
    with open(input_file, 'r') as input_file, open(output_file, 'w') as output_file:
        for line in input_file:
            fields = line.strip().split('\t')
            first_field = fields[0]
            comma_index = first_field.find(',')
            new_first_field = first_field[:comma_index]
            if len(new_first_field) > 13:
                new_first_field = new_first_field[-13:]
            tab_separated_fields = "\t".join(fields[1:])
            output_line = f"{new_first_field}\t{tab_separated_fields}"
            output_file.write(output_line + '\n')

description_process('MC743.fna', 'MC743.fsa')
genemark_process('MC743.genemark.gtf', 'MC743.gtf')
