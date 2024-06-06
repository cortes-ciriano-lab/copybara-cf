# COPYBARA
Copy number tool for ONT cfDNA

## Development / To Do

### Priorities
1. Multiprocessing by chromosome (read counting, CBS and merging) to speed up and avoid merging over chromosomes 
2. Add in option to skip multiprocessing if threads = 1
3. add arg for different normalisation options [self | pon | norm]

### Later
1. fix blacklisting
2. fix centromere regions
3. fix ways chromosomes are defined (e.g. separate ref/resource file)
4. optimise parameters
5. improve plotting function to display fewer points for smaller bin sizes

### Downstream implementations for cfDNA tool
1. tMAD
2. Rascal CN
3. pull out CN features

### Testing/validation
1. COLO829 simulations (dilution series split reads)
2. ONT savana cohort high purity simulations (dilution series split reads)
3. Comparison of ONT COPYBARA vs Illumina ichorCNA


