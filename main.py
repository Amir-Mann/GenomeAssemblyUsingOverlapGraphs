import random
import argparse
from strands_graph import Graph


def read_fasta(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    return ''.join(line.strip() for line in lines if not line.startswith('>'))


def generate_reads(genome, num_reads, read_length, error_prob):
    genome_length = len(genome)
    bases = ['A', 'C', 'G', 'T']
    
    reads = []
    for _ in range(num_reads):
        start_idx = random.randint(0, genome_length - read_length)
        read = list(genome[start_idx:start_idx + read_length])
        
        if error_prob > 0:
            #print("here")
            for j in range(read_length):
                if random.random() < error_prob:
                    read[j] = random.choice([b for b in bases if b != read[j]])
        
        reads.append("".join(read))
    
    return reads


def main():
    parser = argparse.ArgumentParser(description="Generate error-prone reads from a genome.")
    parser.add_argument("-f", "--fasta", type=str, default="demo.fasta", help="Path to the FASTA file containing the genome sequence.")
    parser.add_argument("-l", "--read_length", type=int, default=100, help="Length of each read (default: 100).")
    parser.add_argument("-c", "--coverage", type=float, default=50.0, help="Desired coverage (default: 50x).")
    parser.add_argument("-p", "--error_prob", type=float, default=0.01, help="Base mismatch probability (default: 0.01).")
    args = parser.parse_args()
    
    genome = read_fasta(args.fasta)
    genome_length = len(genome)
    
    num_reads = int((args.coverage * genome_length) / args.read_length)
    
    clean_reads = generate_reads(genome, num_reads, args.read_length, 0)
    
    error_reads = generate_reads(genome, num_reads, args.read_length, args.error_prob)
    
    G = Graph()
    G.load_from_strands(clean_reads)
    sequence = G.get_sequenced_result()
    i = 0
    step = 20
    while i + step < len(genome) and i + step < len(sequence):
        print(f"org: {genome[i: i+step]}")
        print(f"seq: {sequence[i: i+step]}, {sequence[i: i+step]==genome[i: i+step]}")
        i = i + step
    print(f"org: {genome[i:]}")
    print(f"seq: {sequence[i:]}")

    """print(f"Generated {num_reads} reads for {args.coverage}x coverage with error probability {args.error_prob}.")
    print(f"Min clean read is {min(clean_reads)}")
    print(f"Min error read is {min(error_reads)}")"""
    

if __name__ == "__main__":
    main()
