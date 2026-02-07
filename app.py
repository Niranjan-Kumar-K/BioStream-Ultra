from flask import Flask, render_template, request
import re

app = Flask(__name__)

def analyze_dna(seq):
    # Clean the sequence
    seq = seq.upper().strip().replace("\r", "").replace("\n", "").replace(" ", "")
    
    # Validation
    if not all(base in "ATGCN" for base in seq) or len(seq) == 0:
        return "error"

    length = len(seq)
    counts = {base: seq.count(base) for base in "ATCG"}
    
    # 1. GC Content
    gc_val = (seq.count('G') + seq.count('C')) / length * 100

    # 2. Reverse Complement
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    rev_comp = "".join(complement.get(base, base) for base in reversed(seq))

    # 3. Molecular Weight (Approximate)
    mw_dna = (counts['A']*313.2) + (counts['T']*304.2) + (counts['C']*289.2) + (counts['G']*329.2) + 79.0
    
    # 4. Translation Logic
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }
    protein = ""
    for i in range(0, len(seq) - (len(seq) % 3), 3):
        protein += codon_table.get(seq[i:i+3], '?')
    
    mw_protein = len(protein.replace("_","")) * 0.11 # in kDa

    # 5. Restriction Mapping
    enzymes = {"EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT"}
    cuts = {}
    for name, site in enzymes.items():
        positions = [m.start() + 1 for m in re.finditer(f'(?={site})', seq)]
        cuts[name] = positions if positions else ["None"]

    # 6. ORF Finder (Simple ATG search)
    orfs = [m.start() + 1 for m in re.finditer('ATG', seq)]

    return {
        "length": length,
        "counts": counts,
        "gc_content": f"{gc_val:.2f}%",
        "mw_dna": f"{mw_dna/1000:.2f} kDa",
        "mw_protein": f"{mw_protein:.2f} kDa",
        "rev_comp": rev_comp,
        "protein": protein,
        "cuts": cuts,
        "orfs": orfs[:10] # Top 10 ORFs
    }

@app.route("/", methods=["GET", "POST"])
def index():
    results = None
    if request.method == "POST":
        sequence = request.form.get("sequence", "")
        results = analyze_dna(sequence)
    return render_template("index.html", results=results)

if __name__ == "__main__":
    app.run(port=5050, debug=True)