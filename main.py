"""
Main module for assembling the PhiX genome.

This script generates reads from a genome, builds an overlap graph using error-prone or clean reads,
and evaluates the assembly performance. It uses the Graph and timer functions from strands_graph.py.
"""

import math
import time
import tqdm
import random
import argparse
from strands_graph import Graph, timings, timer


def read_fasta(filename: str) -> str:
    """
    Reads a FASTA file and returns the genome sequence as a single string.
    """
    with open(filename, 'r') as file:
        lines = file.readlines()
    return ''.join(line.strip() for line in lines if not line.startswith('>'))


def generate_reads(genome: str, num_reads: int, read_length: int, error_prob: float) -> list[str]:
    """
    Generates a list of reads from the genome sequence, optionally introducing errors.
    """
    genome_length = len(genome)
    bases = ['A', 'C', 'G', 'T']
    
    reads = []
    for _ in range(num_reads):
        start_idx = random.randint(0, genome_length - read_length)
        read = list(genome[start_idx:start_idx + read_length])
        
        if error_prob > 0:
            for j in range(read_length):
                if random.random() < error_prob:
                    read[j] = random.choice([b for b in bases if b != read[j]])
        
        reads.append("".join(read))
    
    return reads


def global_alignment(genome: str, sequence: str) -> tuple[float, float, float]:
    """
    Computes a global alignment between the genome and a given sequence.
    Includes an early stopping mechanism if the alignment is sufficiently bad (e.g., when the offset is 5
    and we are looking at offset 4, most characters would be misaligned).
    """
    genome_length = len(genome)
    log_genome_length = math.log(genome_length, 2)
    
    for i in range(genome_length):
        count = 0
        for j, c in enumerate(sequence):
            if i + j >= genome_length:
                break
            count += int(c == genome[i + j])
            if count - j > log_genome_length * 2:  # Early stopping if alignment is sufficiently bad
                break
        if count > len(sequence) / 3: # If alignment is sufficiently good, use it as the global alignment
            recall = float(count) / genome_length
            precision = float(count) / len(sequence)
            iou = float(count) / (len(sequence) + genome_length - count)
            return recall, precision, iou
    return 0., 0., 0. # Failed to find any good alignments


def iterative_test_sequencing(args: argparse.Namespace, experiment_title: str, genome: str, num_reads: int,
                              allow_mis_matches: int, testing_iterations: int = 1, error_prob: float = 0.0,
                              print_results: bool = True) -> None:
    """
    Iteratively tests the genome assembly process by generating reads, assembling them into a graph,
    and evaluating performance.
    """
    if not args.hide_progress_bar:
        progress_bar = tqdm.tqdm(
            total=testing_iterations,
            desc="Processing",
            bar_format="{l_bar}{bar} {n_fmt}/{total_fmt} | {elapsed} elapsed, {remaining} remaining, {rate_fmt}"
        )

    recall = 0
    precision = 0
    iou = 0
    for _ in range(testing_iterations):
        if not args.hide_progress_bar:
            progress_bar.update(1)
        error_reads = generate_reads(genome, num_reads, args.read_length, error_prob)
        G = Graph()
        with timer("create graph"):
            G.load_from_strands(error_reads, allow_mis_matches)
        
        with timer("sequence"):
            sequence = G.get_sequenced_result()
        
        with timer("eval"):
            r, p, i = global_alignment(genome, sequence)
        
        recall += r / testing_iterations
        precision += p / testing_iterations
        iou += i / testing_iterations

    if not args.hide_progress_bar:
        progress_bar.close()
    
    if print_results:
        print()
        print(experiment_title, "results:")
        print(f"recall: {recall}")
        print(f"precision: {precision}")
        print(f"iou: {iou}")


def main() -> None:
    """
    Main function to parse arguments, generate reads, assemble the genome, and evaluate performance.
    """
    parser = argparse.ArgumentParser(description="Generate error-prone reads from a genome.")
    parser.add_argument("-f", "--fasta", type=str, default="sequence.fasta", help="Path to the FASTA file containing the genome sequence.")
    parser.add_argument("-l", "--read_length", type=int, default=100, help="Length of each read (default: 100).")
    parser.add_argument("-N", "--num_reads", type=int, default=3000, help="Desired amount of reads.")
    parser.add_argument("-p", "--error_prob", type=float, default=0.01, help="Base mismatch probability (default: 0.01).")
    parser.add_argument("-i", "--testing_iterations", type=int, default=100, help="Number of iterations to test the quality of the algorithm. Only for erroneous reads")
    parser.add_argument("-a", "--allow_mis_matches", type=str, choices=["0", "1", "2", "3", "4", "log"], default="2", 
                        help="How many miss matches to allow when comparing suffix to prefix. Only for erroneous reads")

    parser.add_argument("--hide_progress_bar", action="store_true", help="Hides the progress bar")
    parser.add_argument("--hide_timing", action="store_true", help="Hides the timing resultsprogress bar")
    parser.add_argument("-s", "--seed", type=int, default=207732132, help="The seed to run the program with")
    args = parser.parse_args()

    random.seed(args.seed)
    
    genome = read_fasta(args.fasta)
    
    if args.allow_mis_matches == "log":
        allow_mis_matches = int(math.log(args.read_length, 4))
    else:
        allow_mis_matches = int(args.allow_mis_matches)
    
    start = time.time()
    
    title = f"Running with -l {args.read_length} -N {args.num_reads} -p {args.error_prob} -a {args.allow_mis_matches}"
    iterative_test_sequencing(args, title, genome, args.num_reads, allow_mis_matches=allow_mis_matches,
                              testing_iterations=args.testing_iterations, error_prob=args.error_prob)

    total = time.time() - start
    if not args.hide_timing:
        print()
        print(f"took {total:.2f} seconds to run")
        for key, value in timings.items():
            print(f"{key} was {value / total * 100:.2f}% of runtime, taking {value / args.testing_iterations:.3f} seconds to run on average.")
    

if __name__ == "__main__":
    main()
