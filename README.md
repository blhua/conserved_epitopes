# conserved_epitopes
A python script that takes in a Clustal multiple sequence alignment and extracts MHC Class I length epitopes with high conservation across the sequences.

This tool uses a Clustal multiple sequence alignment for any protein of interest, extracts 8-11mer tiling peptides across the alignment, and computes the conservation of the epitope across all the sequences in the alignment. A histogram of the number of epitopes vs. conservation is return along with a csv of all epitopes with a conservation of greater than or equal to 90%.

![conservation](https://github.com/blhua/conserved_epitopes/assets/66856632/e81891b6-6d84-4d90-a188-9720333dfbd8)
![F csv_conservation_plot](https://github.com/blhua/conserved_epitopes/assets/66856632/27938094-0e04-4a88-85e7-e78b6967093e)

# Usage
```python3 get_conserved_epitopes.py [PATH TO ALIGNMENT FILE]```

# Inputs
- Path to alignment file: Path to a Clustal multiple sequence alignment file

# Outputs
- **conservation_plot.jpg**: A jpg file of the histogram of epitope count versus conservation
- **conserved_epitopes.csv**: A comma-separated value file of all epitopes with conservations greater than or equal to 90%
