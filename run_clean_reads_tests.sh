echo "Running tests, this should take about 10-30 min"
echo ""
echo "### Reference Runs ###"
echo "-p 0.0 -a 0 -l 100 -N 3000 -i 100"
python3 main.py -p 0.0 -a 0 -l 100 -N 3000 -i 100
echo "-p 0.0 -a 2 -l 100 -N 3000 -i 100"
python3 main.py -p 0.0 -a 2 -l 100 -N 3000 -i 100
echo ""
echo "### Effect of error probability (p) ###"
echo "-p 0.0 -a 1 -l 100 -N 3000 -i 100"
python3 main.py -p 0.0 -a 1 -l 100 -N 3000 -i 100
echo "-p 0.0 -a 2 -l 100 -N 3000 -i 100"
python3 main.py -p 0.0 -a 2 -l 100 -N 3000 -i 100
echo ""
echo "### Effect of Number of Reads (N) ###"
echo "-p 0.0 -a 2 -l 100 -N 2000 -i 100"
python3 main.py -p 0.0 -a 2 -l 100 -N 2000 -i 100
echo "-p 0.0 -a 2 -l 100 -N 4000 -i 100"
python3 main.py -p 0.0 -a 2 -l 100 -N 4000 -i 100
echo ""
echo "### Effect of Read Length (l) ###"
echo "-p 0.0 -a 2 -l 50 -N 3000 -i 100"
python3 main.py -p 0.0 -a 2 -l 50 -N 3000 -i 100
echo "-p 0.0 -a 2 -l 150 -N 3000 -i 100"
python3 main.py -p 0.0 -a 2 -l 150 -N 3000 -i 100
echo ""
echo "### Effect of number of allowed misalignments (k) ###"
echo "-p 0.0 -a 0 -l 100 -N 3000 -i 100"
python3 main.py -p 0.0 -a 0 -l 100 -N 3000 -i 100
echo "-p 0.0 -a 1 -l 100 -N 3000 -i 100"
python3 main.py -p 0.0 -a 1 -l 100 -N 3000 -i 100
echo "-p 0.0 -a 3 -l 100 -N 3000 -i 100"
python3 main.py -p 0.0 -a 3 -l 100 -N 3000 -i 100