import csv

def process_edges(input_csv, output_txt):
    edges = set()  # Using a set to avoid duplicates
    
    # Read edges from CSV (trying both comma and tab delimiters)
    with open(input_csv, 'r') as csvfile:
        # Try comma delimiter first
        try:
            reader = csv.reader(csvfile, delimiter=',')
            next(reader)  # Skip header if present
            for row in reader:
                if len(row) < 2:
                    continue
                u, v = int(row[0].strip()), int(row[1].strip())
                edges.add((min(u, v), max(u, v)))
        except ValueError:
            # If comma fails, try tab delimiter
            csvfile.seek(0)  # Rewind file
            reader = csv.reader(csvfile, delimiter='\t')
            next(reader)
            for row in reader:
                u, v = int(row[0].strip()), int(row[1].strip())
                edges.add((min(u, v), max(u, v)))
    
    # Write cleaned edges to text file
    with open(output_txt, 'w') as outfile:
        for u, v in sorted(edges):
            outfile.write(f"{u} {v}\n")

if __name__ == "__main__":
    import os 
    
    input_csv = os.path.join("..", "dataset", "politicians", "edges.csv")
    output_txt = os.path.join("..", "dataset", "politicians", "politicians.txt")

    print("Processing file...")
    process_edges(input_csv, output_txt)
    print(f"Successfully processed edges written to {output_txt}")