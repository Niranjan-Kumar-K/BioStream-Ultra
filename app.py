from flask import Flask, render_template, request
from Bio.Seq import Seq
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Data import IUPACData

app = Flask(__name__)

def get_3_letter_protein(one_letter_seq):
    if not one_letter_seq: return "None"
    return "-".join([IUPACData.protein_letters_1to3.get(aa, aa).capitalize() for aa in one_letter_seq])

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/analyzer', methods=['GET', 'POST'])
def analyzer():
    results = None
    if request.method == 'POST':
        seq = request.form.get('sequence', '').strip().upper()
        if seq:
            seq_obj = Seq(seq)
            protein_seq = seq_obj.translate(to_stop=True)
            prot_mass = "0.00 Da"
            if len(protein_seq) > 0:
                prot_mass = f"{ProteinAnalysis(str(protein_seq)).molecular_weight():.2f} Da"
            results = {
                "length": len(seq_obj),
                "dna_mass": f"{molecular_weight(seq_obj, 'DNA'):.2f} Da",
                "protein_mass": prot_mass,
                "translation": get_3_letter_protein(str(protein_seq))
            }
    return render_template('analyzer.html', results=results)

@app.route('/amr', methods=['GET', 'POST'])
def amr():
    amr_results = None
    if request.method == 'POST':
        seq = request.form.get('sequence', '').strip().upper()
        db = {"TGGTATGTGGAAGTTAGATTG": "mecA (Methicillin Resistance)", "TTCGGCATTTCGTC": "vanA (Vancomycin Resistance)"}
        amr_results = next((desc for marker, desc in db.items() if marker in seq), "No known markers found.")
    return render_template('amr.html', amr_results=amr_results)

@app.route('/primer', methods=['GET', 'POST'])
def primer():
    primer_results = None
    if request.method == 'POST':
        seq = request.form.get('sequence', '').strip().upper()
        if len(seq) >= 20:
            f_seq = seq[:20]
            r_seq = str(Seq(seq[-20:]).reverse_complement())
            primer_results = {
                "f_primer": f_seq, "f_tm": f"{molecular_weight(Seq(f_seq), 'DNA')/60:.1f}",
                "r_primer": r_seq, "r_tm": f"{molecular_weight(Seq(r_seq), 'DNA')/60:.1f}"
            }
    return render_template('primer.html', primer_results=primer_results)

if __name__ == '__main__':
    app.run(debug=True)