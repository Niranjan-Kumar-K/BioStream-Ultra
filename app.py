import os, io
from flask import Flask, render_template, request
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis

app = Flask(__name__)

def get_clean_sequence(req):
    if 'fasta_file' in req.files and req.files['fasta_file'].filename != '':
        try:
            file = req.files['fasta_file']
            fasta_data = file.read().decode("utf-8")
            record = SeqIO.read(io.StringIO(fasta_data), "fasta")
            return str(record.seq).upper().strip()
        except: return ""
    raw_seq = req.form.get("sequence", "").upper()
    return "".join(raw_seq.split())

@app.route("/")
def home(): return render_template("index.html")

@app.route("/analyzer", methods=["GET", "POST"])
def analyzer():
    results = None
    if request.method == "POST":
        seq_str = get_clean_sequence(request)
        if seq_str:
            my_seq = Seq(seq_str.replace("U", "T"))
            translation = my_seq.translate(to_stop=True)
            clean_prot = "".join([a for a in str(translation) if a in "ACDEFGHIKLMNPQRSTVWY"])
            p_mass = round(ProteinAnalysis(clean_prot).molecular_weight(), 2) if clean_prot else 0
            
            enzymes = {"GAATTC": "EcoRI", "GGATCC": "BamHI", "AAGCTT": "HindIII", "GGTACC": "KpnI"}
            sites = [f"{n} ({s}) at pos {str(my_seq).find(s)+1}" for s, n in enzymes.items() if s in str(my_seq)]

            results = {
                "length": len(my_seq), "gc_content": round(gc_fraction(my_seq)*100, 2),
                "dna_mass": f"{len(my_seq)*660:,} Da", "prot_mass": f"{p_mass:,} Da",
                "translation": str(translation), "reverse_seq": str(my_seq.reverse_complement()),
                "orfs": [i+1 for i in range(len(my_seq)-2) if my_seq[i:i+3] == "ATG"],
                "restriction_sites": sites if sites else ["None Detected"]
            }
    return render_template("analyzer.html", results=results)

@app.route("/amr", methods=["GET", "POST"])
def amr_scout():
    amr_results = None
    if request.method == "POST":
        seq_str = get_clean_sequence(request)
        db = {"mecA (Methicillin Resistance)": "TGGTATGTGGAAGTTAGATTG", "blaCTX-M": "ATGGTGACAAAGAGAGTGCAA"}
        found = [g for g, s in db.items() if s in seq_str]
        amr_results = "⚠️ WARNING: " + ", ".join(found) if found else "✅ No markers detected."
    return render_template("amr.html", amr_results=amr_results)

@app.route("/primer", methods=["GET", "POST"])
def primer_lab():
    res = None
    if request.method == "POST":
        seq = get_clean_sequence(request)
        if len(seq) >= 50:
            f, r = seq[:20], str(Seq(seq[-20:]).reverse_complement())
            tm = lambda p: round(64.9 + 41 * (p.count('G') + p.count('C') - 16.4) / 20, 1)
            res = {"f_primer": f, "r_primer": r, "f_tm": tm(f), "r_tm": tm(r)}
        else: res = "Sequence too short (min 50bp)."
    return render_template("primer.html", primer_results=res)

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=int(os.environ.get("PORT", 5000)))