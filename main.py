import math
import time
import tqdm
import random
import argparse
from strands_graph import Graph, timings, timer


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
    parser.add_argument("-i", "--testing_iterations", type=int, default=1000, help="Number of iterations to test the quality of the algorithm.")
    parser.add_argument("-a", "--allow_mis_matches", type=str, choices=["0", "1", "2", "3", "log"], default=1000, help="How many miss matches to allow when comparing suffix to prefix.")

    parser.add_argument("--hide_progress_bar", action="store_true", help="Hides the progress bar")
    parser.add_argument("--hide_timing", action="store_true", help="Hides the timing resultsprogress bar")
    parser.add_argument("-s", "--seed", type=int, default=7, help="The seed to run the program with")
    args = parser.parse_args()

    random.seed(args.seed)
    
    genome = read_fasta(args.fasta)
    genome_length = len(genome)
    log_genome_length = math.log(genome_length, 2)
    
    num_reads = int((args.coverage * genome_length) / args.read_length)
    print(f"Generating {num_reads} reads.")
    
    clean_reads = generate_reads(genome, num_reads, args.read_length, 0)

    if args.allow_mis_matches == "log":
        allow_mis_matches = int(math.log(args.read_length, 4))
    else:
        allow_mis_matches = int(args.allow_mis_matches)
    
    G = Graph()
    G.load_from_strands(clean_reads, allow_mis_matches)
    sequence = G.get_sequenced_result()
    for i in range(30):
        if sequence == genome[i:i+len(sequence)]:
            print("match of length", len(sequence), "out of", len(genome))
            break
    else:
        print("No match")
    

    if not args.hide_progress_bar:
        progress_bar = tqdm.tqdm(
            total=args.testing_iterations,
            desc="Processing",
            bar_format="{l_bar}{bar} {n_fmt}/{total_fmt} | {elapsed} elapsed, {remaining} remaining, {rate_fmt}"
        )
    start = time.time()
    
    recall = 0
    precision = 0
    iou = 0
    for _ in range(args.testing_iterations):
        if not args.hide_progress_bar:
            progress_bar.update(1)
        error_reads = generate_reads(genome, num_reads, args.read_length, args.error_prob)
        G = Graph()
        with timer("create graph"):
            G.load_from_strands(error_reads, allow_mis_matches)
        
        with timer("sequence"):
            sequence = G.get_sequenced_result()
        
        with timer("eval"):
            for i in range(genome_length):
                count = 0
                for j, c in enumerate(sequence):
                    if i + j >= genome_length:
                        break
                    count += int(c == genome[i + j])
                    if count - j > log_genome_length * 2: # Early stopping, see report
                        break
                if count > len(sequence) / 3:
                    #print("match of length", len(sequence), "out of", genome_length, "with count", count)
                    recall += float(count) / args.testing_iterations / genome_length
                    precision += float(count) / args.testing_iterations / len(sequence)
                    iou += float(count) / (len(sequence) + genome_length - count) / args.testing_iterations
                    break
    if not args.hide_progress_bar:
        progress_bar.close()
    
    print()
    print(f"recall: {recall}")
    print(f"precision: {precision}")
    print(f"iou: {iou}")
    
    total = time.time() - start
    if not args.hide_timing:
        print(f"took {total:.2f} seconds to run")
        for key, value in timings.items():
            print(f"{key} was {value / total * 100:.2f}% of runtime")
    

if __name__ == "__main__":
    main()
