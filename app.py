import os
import io
from flask import Flask, render_template, request
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

app = Flask(__name__)

# --- HELPER FUNCTIONS ---

def get_clean_sequence(req):
    """Extracts sequence from either file upload or text area."""
    if 'fasta_file' in req.files and req.files['fasta_file'].filename != '':
        file = req.files['fasta_file']
        fasta_data = file.read().decode("utf-8")
        fasta_io = io.StringIO(fasta_data)
        try:
            record = SeqIO.read(fasta_io, "fasta")
            return str(record.seq).upper()
        except:
            return "" # Handle non-fasta files gracefully
    else:
        raw_seq = req.form.get("sequence", "").upper()
        return "".join(raw_seq.split())

# --- ROUTES ---

@app.route("/")
def home():
    return render_template("index.html")

@app.route("/analyzer", methods=["GET", "POST"])
def analyzer():
    results = None
    if request.method == "POST":
        seq_str = get_clean_sequence(request)
        clean_seq = "".join([c for c in seq_str if c in "ATCGU"])
        
        if clean_seq:
            my_seq = Seq(clean_seq.replace("U", "T"))
            results = {
                "length": len(my_seq),
                "gc_content": round(gc_fraction(my_seq) * 100, 2),
                "translation": str(my_seq.translate(to_stop=True)),
                "orfs": [i for i in range(len(clean_seq)-2) if clean_seq[i:i+3] == "ATG"]
            }
    return render_template("analyzer.html", results=results)

@app.route("/amr", methods=["GET", "POST"])
def amr_scout():
    amr_results = None
    if request.method == "POST":
        seq_str = get_clean_sequence(request)
        
        # Starter Dictionary of AMR Gene Signatures (Short motifs for demo)
        # In a real app, these would be 1000+ base pairs long
        amr_database = {
            "mecA (Methicillin Resistance)": "TGGTATGTGGAAGTTAGATTG",
            "blaCTX-M (Beta-lactamase)": "ATGGTGACAAAGAGAGTGCAA",
            "tetA (Tetracycline Resistance)": "GTGAAACCCAACATACCCCCC",
            "vanA (Vancomycin Resistance)": "CATGAATAGAATAAAAGTTGC"
        }
        
        found_genes = []
        for gene_name, signature in amr_database.items():
            if signature in seq_str:
                found_genes.append(gene_name)
        
        if found_genes:
            amr_results = "⚠️ WARNING: Matches found for: " + ", ".join(found_genes)
        else:
            amr_results = "✅ No high-priority AMR markers detected in the provided sequence."
            
    return render_template("amr.html", amr_results=amr_results)

@app.route("/primer", methods=["GET", "POST"])
def primer_lab():
    primer_results = None
    if request.method == "POST":
        seq_str = get_clean_sequence(request)
        
        if len(seq_str) > 50:
            # Simple Primer Design: 20bp from start and 20bp from end
            f_primer = seq_str[:20]
            # Reverse primer must be the reverse complement of the end
            r_primer_raw = Seq(seq_str[-20:])
            r_primer = str(r_primer_raw.reverse_complement())
            
            primer_results = {
                "f_primer": f_primer,
                "r_primer": r_primer,
                "f_tm": round(64.9 + 41 * (f_primer.count('G') + f_primer.count('C') - 16.4) / len(f_primer), 1),
                "r_tm": round(64.9 + 41 * (r_primer.count('G') + r_primer.count('C') - 16.4) / len(r_primer), 1)
            }
        else:
            primer_results = "Error: Sequence too short for primer design (min 50bp)."

    return render_template("primer.html", primer_results=primer_results)

if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)