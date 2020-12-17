#!/usr/bin/env py

import sys

sites_vcf = sys.argv[1]

alleles = {
'0/0' : 'ref',
'1/1' : 'alt',
'./.' : 'N',
'.|.' : 'N',
'0/1' : '.',
'0|1' : '.',
'1|0' : '.',
'1|1' : 'alt',
'0|0' : 'ref'
}

previous_chromosome = 'BEGIN_FILE'
output = ''
with open(sites_vcf, 'r') as infile:
    for line in infile:
        if not line.startswith("##"):
            splitline = line.rstrip().split()
            if line.startswith("#"):
                line_ids = splitline[9:]
                n_lines = len(line_ids)
                continue
            else:
                chromosome = splitline[0]
                position = splitline[1]
                alleles['0/0'] = splitline[3]
                alleles['0|0'] = splitline[3]
                alleles['1/1'] = splitline[4]
                alleles['1|1'] = splitline[4]
                genotypes = [alleles[x] for x in splitline[9:]]
            if previous_chromosome == 'BEGIN_FILE' :
                previous_chromosome = chromosome
                outfile = open('''%s.H12.csv''' % (chromosome), 'w')
                print("Working on %s.H12.csv" % (chromosome))
                #header = ','.join([str(x) for x in ([chromosome, 'Ref'] + line_ids + ["Coverage"])]) + '\n'
                #outfile.write(header)
            if previous_chromosome == chromosome:
                vals = [position] + genotypes
                row = ','.join([str(x) for x in vals]) + ',\n'
                outfile.write(row)
            else:
                outfile.close()
                previous_chromosome = chromosome
                outfile = open('''%s.H12.csv''' % (chromosome), 'w')
                print("Working on %s.H12.csv" % (chromosome))
                #header = ','.join([str(x) for x in ([chromosome, 'Ref'] + line_ids + ["Coverage"])]) + '\n'
                #outfile.write(header)
    outfile.close()
