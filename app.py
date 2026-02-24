import os
import io
from flask import Flask, render_template, request
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.SeqUtils.ProtParam import ProteinAnalysis

app = Flask(__name__)

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
            translation = my_seq.translate(to_stop=True)
            
            # --- CALCULATIONS ---
            # DNA Mass: Approx 660 Daltons per base pair
            dna_mass = len(my_seq) * 660 
            
            # Protein Mass
            prot_mass = 0
            if len(translation) > 0:
                # Remove any characters that aren't standard amino acids for the calculator
                clean_prot = "".join([amino for amino in str(translation) if amino in "ACDEFGHIKLMNPQRSTVWY"])
                if clean_prot:
                    analysed_seq = ProteinAnalysis(clean_prot)
                    prot_mass = round(analysed_seq.molecular_weight(), 2)

            # --- RESTRICTION MAPPING ---
            enzymes = {"GAATTC": "EcoRI", "GGATCC": "BamHI", "AAGCTT": "HindIII", "GGTACC": "KpnI"}
            found_sites = []
            for site, name in enzymes.items():
                pos = str(my_seq).find(site)
                while pos != -1:
                    found_sites.append(f"{name} ({site}) at pos {pos + 1}")
                    pos = str(my_seq).find(site, pos + 1)

            results = {
                "length": len(my_seq),
                "gc_content": round(gc_fraction(my_seq) * 100, 2),
                "dna_mass": f"{dna_mass:,} Da",
                "prot_mass": f"{prot_mass:,} Da",
                "translation": str(translation),
                "reverse_seq": str(my_seq.reverse_complement()),
                "orfs": [i + 1 for i in range(len(my_seq)-2) if my_seq[i:i+3] == "ATG"],
                "restriction_sites": found_sites if found_sites else ["No common sites found"]
            }
    return render_template("analyzer.html", results=results)

# (Keep your /amr and /primer routes as they were!)

if __name__ == "__main__":
    app.run(host='0.0.0.0', port=int(os.environ.get("PORT", 5000)))