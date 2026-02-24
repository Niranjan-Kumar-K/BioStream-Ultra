import os
import io
from flask import Flask, render_template, request
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction

app = Flask(__name__)

# --- HELPER FUNCTIONS ---
def get_clean_sequence(req):
    if 'fasta_file' in req.files and req.files['fasta_file'].filename != '':
        file = req.files['fasta_file']
        fasta_data = file.read().decode("utf-8")
        fasta_io = io.StringIO(fasta_data)
        try:
            record = SeqIO.read(fasta_io, "fasta")
            return str(record.seq).upper()
        except:
            return ""
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
            # Handle DNA and RNA
            my_seq = Seq(clean_seq.replace("U", "T"))
            
            # Restriction Enzyme Logic
            enzymes = {
                "GAATTC": "EcoRI",
                "GGATCC": "BamHI",
                "AAGCTT": "HindIII",
                "GGTACC": "KpnI"
            }
            
            found_sites = []
            for site, name in enzymes.items():
                pos = str(my_seq).find(site)
                while pos != -1:
                    found_sites.append(f"{name} ({site}) at pos {pos + 1}")
                    pos = str(my_seq).find(site, pos + 1)

            results = {
                "length": len(my_seq),
                "gc_content": round(gc_fraction(my_seq) * 100, 2),
                "translation": str(my_seq.translate(to_stop=True)),
                "orfs": [i + 1 for i in range(len(my_seq)-2) if my_seq[i:i+3] == "ATG"],
                "restriction_sites": found_sites if found_sites else ["No common sites found"]
            }
    return render_template("analyzer.html", results=results)

@app.route("/amr", methods=["GET", "POST"])
def amr_scout():
    amr_results = None
    if request.method == "POST":
        seq_str = get_clean_sequence(request)
        amr_database = {
            "mecA (Methicillin Resistance)": "TGGTATGTGGAAGTTAGATTG",
            "blaCTX-M (Beta-lactamase)": "ATGGTGACAAAGAGAGTGCAA",
            "tetA (Tetracycline Resistance)": "GTGAAACCCAACATACCCCCC"
        }
        found_genes = [gene for gene, sig in amr_database.items() if sig in seq_str]
        amr_results = "⚠️ WARNING: Matches for: " + ", ".join(found_genes) if found_genes else "✅ No high-priority AMR markers detected."
    return render_template("amr.html", amr_results=amr_results)

@app.route("/primer", methods=["GET", "POST"])
def primer_lab():
    primer_results = None
    if request.method == "POST":
        seq_str = get_clean_sequence(request)
        if len(seq_str) >= 50:
            f_primer = seq_str[:20]
            r_primer = str(Seq(seq_str[-20:]).reverse_complement())
            primer_results = {
                "f_primer": f_primer,
                "r_primer": r_primer,
                "f_tm": round(64.9 + 41 * (f_primer.count('G') + f_primer.count('C') - 16.4) / 20, 1),
                "r_tm": round(64.9 + 41 * (r_primer.count('G') + r_primer.count('C') - 16.4) / 20, 1)
            }
        else:
            primer_results = "Error: Sequence too short (min 50bp)."
    return render_template("primer.html", primer_results=primer_results)

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=int(os.environ.get("PORT", 5000)))