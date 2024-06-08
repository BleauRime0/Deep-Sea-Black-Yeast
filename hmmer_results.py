import csv

def filter_hmmer_output(input_file, output_file, evalue_threshold=1e-10, score_threshold=0):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        tsv_writer = csv.writer(outfile, delimiter='\t')
        tsv_writer.writerow(['Gene Name', 'Accession', 'Protein Sequence ID', 'E value', 'Score', 'Bias', 'Protein Function'])

        for line in infile:
            if line.startswith('#') or not line.strip():
                continue

            parts = line.split()
            target_name = parts[0]
            accession = parts[1]
            query_name = parts[2]
            evalue = float(parts[4])
            score = float(parts[5])
            bias = float(parts[6])
            description = ' '.join(parts[18:])

            if evalue <= evalue_threshold and score-bias >= score_threshold:
                tsv_writer.writerow([target_name, accession, query_name, evalue, score, bias, description])

filter_hmmer_output('MC743_COGs_hmmer.txt', 'MC743_protein_function.tsv')
