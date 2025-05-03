def load_output(file_path):
    result = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            node, dist, parent = parts
            result[int(node)] = (dist, parent)
    return result

def compare_outputs(file1, file2):
    output1 = load_output(file1)
    output2 = load_output(file2)

    if output1.keys() != output2.keys():
        print("Mismatch in node set!")
        #return

    mismatches = []
    for node in output1:
        if output1[node] != output2[node]:
            mismatches.append((node, output1[node], output2[node]))

    if mismatches:
        print(f"Found {len(mismatches)} mismatches:")
        for node, val1, val2 in mismatches:
            print(f"Node {node}: File1 = {val1}, File2 = {val2}")
    else:
        print("✅ Outputs are identical!")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python compare_outputs.py output1.txt sssp_result.txt")
        sys.exit(1)

    compare_outputs("output/output1.txt", "output/sssp_result.txt")