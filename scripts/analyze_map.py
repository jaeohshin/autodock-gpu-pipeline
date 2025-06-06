def read_map_file(filename):
    values = []
    with open(filename, 'r') as f:
        lines = f.readlines()[6:]  # skip header
        for line in lines:
            # split by whitespace and filter out any empty strings
            nums = [float(x) for x in line.strip().split() if x]
            values.extend(nums)
    return values

def main():
    import sys
    import numpy as np

    filename = sys.argv[1]
    values = read_map_file(filename)
    values_np = np.array(values)

    print(f"[{filename}]")
    print(f"  Number of points: {len(values_np)}")
    print(f"  Mean energy     : {np.mean(values_np):.3f} kcal/mol")
    print(f"  Median energy   : {np.median(values_np):.3f}")
    print(f"  Max energy      : {np.max(values_np):.3f}")
    print(f"  # Points > 10k  : {(values_np > 10000).sum()}")
    print(f"  # Points > 1e5  : {(values_np > 1e5).sum()}")

if __name__ == "__main__":
    main()
