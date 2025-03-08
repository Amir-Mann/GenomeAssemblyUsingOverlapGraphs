# Genome Assembly Using Overlap Graphs

Tested with Python 3.10.12.

## Requirements
- Python 3.10.12
- tqdm (install via `pip install tqdm`)  
  Alternatively, use the `--hide_progress_bar` flag to disable progress visualization.

## Usage
- **Main Script:**  
  Run the main script to assemble the PhiX genome:  

`
python3 main.py -f sequence.fasta -l 100 -N 3000 -p 0.01 -a 2 -i 100
`


- **Reported Results:**  
The results reported in the project report could be replicated by running:

`
./run_tests.sh
`

- **Clean Reads Testing:**  
To replicate results on error-free (clean) reads, run:
 
`
./run_clean_reads_tests.sh
`
